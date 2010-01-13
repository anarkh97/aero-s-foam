#include "LocalNetwork.h"

#include "../RemoteStateTask.h"
#include "../ZeroSharedState.h"

#include <algorithm>
#include <iterator>

#include <cassert>

namespace Pita { namespace Hts {

/*struct SameParity {
public:
  explicit SameParity(HalfSliceRank ref) :
    ref_(ref)
  {}

  template <typename T>
  bool operator()(const std::pair<const HalfSliceRank, T> p) const {
    return (p.first.value() - ref_.value() % 2) != 0;
  }

private:
  HalfSliceRank ref_;
}*/

class ZeroReducedCorrection : public NamedTask {
public:
  EXPORT_PTRINTERFACE_TYPES(ZeroReducedCorrection);

  void iterationIs(IterationRank i) {
    ZeroSharedState<Vector>::Ptr zeroCorrection = ZeroSharedState<Vector>::New("Impl Zero Reduced State");
    zeroCorrection->targetIs(target_.ptr());
    zeroCorrection->stateSizeIs(assembler_->reducedBasisSize());

    zeroCorrection->iterationIs(i);
  }

  explicit ZeroReducedCorrection(ReducedSeed * target, const UpdatedSeedAssembler * assembler) :
    NamedTask(String("Zero Reduced Correction ") + target->name()),
    target_(target),
    assembler_(assembler)
  {}

private:
  ReducedSeed::Ptr target_;
  UpdatedSeedAssembler::PtrConst assembler_;
};

class ZeroCorrection : public NamedTask {
public:
  EXPORT_PTRINTERFACE_TYPES(ZeroCorrection);

  void iterationIs(IterationRank i) {
    ZeroSharedState<DynamState>::Ptr zeroCorrection = ZeroSharedState<DynamState>::New("Impl Zero State");
    zeroCorrection->targetIs(target_.ptr());
    zeroCorrection->stateSizeIs(assembler_->vectorSize());

    zeroCorrection->iterationIs(i);
  }

  explicit ZeroCorrection(Seed * target, const UpdatedSeedAssembler * assembler) :
    NamedTask(String("Zero Correction ") + target->name()),
    target_(target),
    assembler_(assembler)
  {}

private:
  Seed::Ptr target_;
  UpdatedSeedAssembler::PtrConst assembler_;
};

LocalNetwork::LocalNetwork(SliceMapping * mapping,
                           CpuRank localCpu,
                           HalfTimeSlice::Manager * sliceMgr,
                           JumpProjector::Manager * jumpProjMgr,
                           ReducedFullTimeSlice::Manager * ftsMgr,
                           UpdatedSeedAssembler::Manager * usaMgr,
                           CorrectionTimeSlice::Manager * ctsMgr,
                           RemoteState::Manager * commMgr,
                           SeedErrorEvaluator::Manager * jumpErrorMgr) :
  mapping_(mapping),
  localCpu_(localCpu),
  sliceMgr_(sliceMgr),
  jumpProjMgr_(jumpProjMgr),
  ftsMgr_(ftsMgr),
  usaMgr_(usaMgr),
  ctsMgr_(ctsMgr),
  commMgr_(commMgr),
  seedMgr_(Seed::Manager::New()),
  reducedSeedMgr_(ReducedSeed::Manager::New()),
  jumpErrorMgr_(jumpErrorMgr)
{
  init();
}

void
LocalNetwork::init() {
  for (SliceMapping::SliceIterator it = mapping_->hostedSlice(localCpu()); it; ++it) {
    HalfSliceRank sliceRank = *it;
    HalfSliceRank nextSliceRank = sliceRank + HalfSliceCount(1);
    HalfSliceRank previousSliceRank = sliceRank - HalfSliceCount(1);
    HalfSliceRank previousFullSliceRank = sliceRank - HalfSliceCount(2);
    HalfSliceRank nextFullSliceRank = sliceRank + HalfSliceCount(2);
    
    CpuRank previousCpu = mapping_->hostCpu(previousSliceRank);
    CpuRank nextCpu(mapping_->hostCpu(nextSliceRank));
    CpuRank previousFullSliceCpu(mapping_->hostCpu(previousFullSliceRank));
    CpuRank nextFullSliceCpu(mapping_->hostCpu(nextFullSliceRank));
    
    // HalfTimeSlices
    HalfTimeSlice::Ptr forwardSlice = sliceMgr_->instanceNew(HalfSliceId(sliceRank, HalfTimeSlice::FORWARD));
    forwardSlice->seedIs(getSeed(SeedId(MAIN_SEED, sliceRank)));
    Seed::Ptr nextLeftPropagatedSeed = getSeed(SeedId(LEFT_SEED, nextSliceRank));
    forwardSlice->propagatedSeedIs(nextLeftPropagatedSeed.ptr());
    halfTimeSlice_.push_back(forwardSlice);

    if (sliceRank > HalfSliceRank(0)) {
      HalfTimeSlice::Ptr backwardSlice = sliceMgr_->instanceNew(HalfSliceId(sliceRank, HalfTimeSlice::BACKWARD));
      backwardSlice->seedIs(getSeed(SeedId(MAIN_SEED, nextSliceRank)));
      backwardSlice->propagatedSeedIs(getSeed(SeedId(RIGHT_SEED, sliceRank)));
      halfTimeSlice_.push_back(backwardSlice);
    }

    // JumpProjectors
    JumpProjector::Ptr jumpProjector = jumpProjMgr_->instanceNew(toString(sliceRank));
    jumpProjector->predictedSeedIs(getSeed(SeedId(RIGHT_SEED, sliceRank)));
    Seed::Ptr leftPropagatedSeed = getSeed(SeedId(LEFT_SEED, sliceRank));
    jumpProjector->actualSeedIs(leftPropagatedSeed.ptr());
    Seed::Ptr fullJump = getSeed(SeedId(SEED_JUMP, sliceRank));
    jumpProjector->seedJumpIs(fullJump.ptr());
    jumpProjector->reducedSeedJumpIs(getReducedSeed(SeedId(SEED_JUMP, sliceRank)));
    jumpProjector_.insert(jumpProjector_.end(), std::make_pair(sliceRank, jumpProjector));

    if (jumpErrorMgr_) {
      SeedErrorEvaluator::Ptr jumpErrorEvaluator = jumpErrorMgr_->instanceNew(fullJump.ptr());
      jumpErrorEvaluator->referenceSeedIs(leftPropagatedSeed.ptr());
    }

    // FullTimeSlices
    TaskMap & fullTimeSlice = ((sliceRank.value()) % 2 == 0) ? evenFullTimeSlice_ : oddFullTimeSlice_;
    ReducedSeed::Ptr correction = getReducedSeed(SeedId(SEED_CORRECTION, sliceRank));
    ReducedSeed::Ptr nextCorrection(NULL);
    if (nextSliceRank < firstInactiveSlice()) {
      ReducedFullTimeSlice::Ptr fullSlice = ftsMgr_->instanceNew(sliceRank);
      fullSlice->jumpIs(getReducedSeed(SeedId(SEED_JUMP, sliceRank)));
      fullSlice->correctionIs(correction.ptr());
      nextCorrection = getReducedSeed(SeedId(SEED_CORRECTION, nextFullSliceRank));
      fullSlice->nextCorrectionIs(nextCorrection.ptr());
      //log() << "Setup FTS head/tail = " << fullSlice->headHalfSlice() << "/" << fullSlice->tailHalfSlice() << "\n";
      fullTimeSlice.insert(std::make_pair(sliceRank, fullSlice));
    }

    // Update
    String forwardAssemblerName = toString(sliceRank); 
    UpdatedSeedAssembler::Ptr forwardAssembler = usaMgr_->instance(forwardAssemblerName);
    if (!forwardAssembler) {
      UpdatedSeedAssembler::Ptr forwardAssembler = usaMgr_->instanceNew(forwardAssemblerName);
      forwardAssembler->correctionIs(getSeed(SeedId(SEED_CORRECTION, sliceRank)));
      forwardAssembler->propagatedSeedIs(getSeed(SeedId(LEFT_SEED, sliceRank)));
      forwardAssembler->correctionComponentsIs(getReducedSeed(SeedId(SEED_CORRECTION, sliceRank)));
      forwardAssembler->updatedSeedIs(getSeed(SeedId(MAIN_SEED, sliceRank)));
      seedAssembler_.insert(seedAssembler_.end(), std::make_pair(sliceRank, forwardAssembler));
    }

    String backwardAssemberName = toString(nextSliceRank); // Performs duplicate work for communication optimization
    UpdatedSeedAssembler::Ptr backwardAssembler = usaMgr_->instanceNew(backwardAssemberName);
    backwardAssembler->correctionIs(getSeed(SeedId(SEED_CORRECTION, nextSliceRank)));
    backwardAssembler->propagatedSeedIs(getSeed(SeedId(LEFT_SEED, nextSliceRank)));
    ReducedSeed::Ptr reducedCorrection = getReducedSeed(SeedId(SEED_CORRECTION, nextSliceRank));
    backwardAssembler->correctionComponentsIs(reducedCorrection.ptr());
    backwardAssembler->updatedSeedIs(getSeed(SeedId(MAIN_SEED, nextSliceRank)));
    seedAssembler_.insert(seedAssembler_.end(), std::make_pair(nextSliceRank, backwardAssembler));
    
    // Communications: Left Propagated Seeds
    if (previousCpu != CpuRank(-1) && previousCpu != localCpu()) {
      RemoteState::SeedWriter::Ptr writer = commMgr_->writerNew(leftPropagatedSeed.ptr(), previousCpu);
      if (writer) {
        String taskName = String("Receive Propagated Seed ") + toString(sliceRank);
        RemoteStateTask::Ptr task = RemoteStateTask::New(taskName, writer.ptr());
        leftSeedSync_.insert(std::make_pair(sliceRank, task));
      }
    }
    if (nextCpu != CpuRank(-1) && nextCpu != localCpu()) {
      RemoteState::SeedReader::Ptr reader = commMgr_->readerNew(nextLeftPropagatedSeed.ptr(), nextCpu);
      if (reader) {
        String taskName = String("Send Propagated Seed ") + toString(nextSliceRank);
        RemoteStateTask::Ptr task = RemoteStateTask::New(taskName, reader.ptr());
        leftSeedSync_.insert(std::make_pair(nextSliceRank, task));
      }
    }

    // Communications: Correction during sequential phase
    if (previousFullSliceCpu != CpuRank(-1) && previousFullSliceCpu != localCpu()) {
      RemoteState::ReducedSeedWriter::Ptr writer = commMgr_->writerNew(correction.ptr(), previousFullSliceCpu);
      if (writer) {
        String taskName = String("Receive Correction ") + toString(sliceRank);
        RemoteStateTask::Ptr task = RemoteStateTask::New(taskName, writer.ptr());
        fullTimeSlice.insert(std::make_pair(previousSliceRank, task));
      }
    }
    if (nextCorrection) {    
      if (nextFullSliceCpu != CpuRank(-1) && nextFullSliceCpu != localCpu()) {
        RemoteState::ReducedSeedReader::Ptr reader = commMgr_->readerNew(nextCorrection.ptr(), nextFullSliceCpu);
        if (reader) {
          String taskName = String("Send Correction ") + toString(nextFullSliceRank);
          RemoteStateTask::Ptr task = RemoteStateTask::New(taskName, reader.ptr());
          fullTimeSlice.insert(std::make_pair(nextSliceRank, task));
        }
      }
    }

    // Communication: Correction as an optimization
    if (previousCpu != CpuRank(-1) && localCpu() != previousCpu && localCpu() != nextCpu) {
      RemoteState::ReducedSeedWriter::Ptr correctionWriter = commMgr_->writerNew(reducedCorrection.ptr(), previousCpu);
      if (correctionWriter) {
        String taskName = String("Receive Correction ") + toString(nextSliceRank);
        RemoteStateTask::Ptr task = RemoteStateTask::New(taskName, correctionWriter.ptr());
        correctionSync_.insert(correctionSync_.end(), std::make_pair(previousSliceRank, task));
      }
    }
    if (nextCorrection) {    
      if (nextCpu != CpuRank(-1) && nextCpu != localCpu() && nextCpu != nextFullSliceCpu) {
        RemoteState::ReducedSeedReader::Ptr correctionReader = commMgr_->readerNew(nextCorrection.ptr(), nextCpu);
        if (correctionReader) {
          String taskName = String("Send Correction ") + toString(nextFullSliceRank);
          RemoteStateTask::Ptr task = RemoteStateTask::New(taskName, correctionReader.ptr());
          correctionSync_.insert(correctionSync_.end(), std::make_pair(sliceRank, task));
        }
      }
    }

    // Affine (HACK)
    if (ctsMgr_) {
      TaskMap & coarseTimeSlice = ((sliceRank.value()) % 2 == 0) ? evenCoarseTimeSlice_ : oddCoarseTimeSlice_;
      
      Seed::Ptr correction = getSeed(SeedId(SEED_CORRECTION, sliceRank));
      Seed::Ptr nextCorrection(NULL);
      // Coarse Propagation
      if (nextSliceRank < firstInactiveSlice()) {
        CorrectionTimeSlice::Ptr correctionSlice = ctsMgr_->instanceNew(sliceRank);
        correctionSlice->predictedSeedIs(getSeed(SeedId(RIGHT_SEED, sliceRank)));
        correctionSlice->actualSeedIs(getSeed(SeedId(LEFT_SEED, sliceRank)));
        correctionSlice->correctionIs(correction.ptr());
        correctionSlice->jumpIs(getSeed(SeedId(SEED_JUMP, sliceRank)));
        nextCorrection = getSeed(SeedId(SEED_CORRECTION, nextFullSliceRank));
        correctionSlice->nextCorrectionIs(nextCorrection.ptr());
        coarseTimeSlice.insert(std::make_pair(sliceRank, correctionSlice.ptr()));
      }

      // Communications: Correction during coarse propagation (sequential phase)
      if (previousFullSliceCpu != CpuRank(-1) && previousFullSliceCpu != localCpu()) {
        RemoteState::SeedWriter::Ptr writer = commMgr_->writerNew(correction.ptr(), previousFullSliceCpu);
        if (writer) {
          String taskName = String("Receive (Affine) Coarse Correction ") + toString(sliceRank);
          RemoteStateTask::Ptr task = RemoteStateTask::New(taskName, writer.ptr());
          coarseTimeSlice.insert(std::make_pair(previousFullSliceRank, task));
        }
      }
      if (nextCorrection) {    
        if (nextFullSliceCpu != CpuRank(-1) && nextFullSliceCpu != localCpu()) {
          RemoteState::SeedReader::Ptr reader = commMgr_->readerNew(nextCorrection.ptr(), nextFullSliceCpu);
          if (reader) {
            String taskName = String("Send (Affine) Coarse Correction ") + toString(nextFullSliceRank);
            RemoteStateTask::Ptr task = RemoteStateTask::New(taskName, reader.ptr());
            coarseTimeSlice.insert(std::make_pair(nextFullSliceRank, task));
          }
        }
      }

      // Communication: Correction as an optimization
      if (previousCpu != CpuRank(-1) && localCpu() != previousCpu && localCpu() != nextCpu) {
        RemoteState::SeedWriter::Ptr correctionWriter = commMgr_->writerNew(getSeed(SeedId(SEED_CORRECTION, nextSliceRank)), previousCpu);
        if (correctionWriter) {
          String taskName = String("Receive (Affine/Optim) Coarse Correction ") + toString(nextSliceRank);
          RemoteStateTask::Ptr task = RemoteStateTask::New(taskName, correctionWriter.ptr());
          fullCorrectionSync_.insert(fullCorrectionSync_.end(), std::make_pair(nextSliceRank, task));
        }
      }
      if (nextCorrection) {    
        if (nextCpu != CpuRank(-1) && nextCpu != localCpu() && nextCpu != nextFullSliceCpu) {
          RemoteState::SeedReader::Ptr correctionReader = commMgr_->readerNew(nextCorrection.ptr(), nextCpu);
          if (correctionReader) {
            String taskName = String("Send (Affine/Optim) Coarse Correction ") + toString(nextFullSliceRank);
            RemoteStateTask::Ptr task = RemoteStateTask::New(taskName, correctionReader.ptr());
            fullCorrectionSync_.insert(fullCorrectionSync_.end(), std::make_pair(nextFullSliceRank, task));
          }
        }
      }
    }
  }
}

Seed *
LocalNetwork::getSeed(const SeedId & id) {
  String seedName = toString(id);
  Seed::Ptr seed = seedMgr_->instance(seedName);
  if (!seed) {
    seed = seedMgr_->instanceNew(seedName);
    if (id.type() == MAIN_SEED) {
      mainSeed_.insert(std::make_pair(id.rank(), seed));
    }
  }

  return seed.ptr();
}

ReducedSeed *
LocalNetwork::getReducedSeed(const SeedId & id) {
  String seedName = toString(id);
  ReducedSeed::Ptr seed = reducedSeedMgr_->instance(seedName);
  if (!seed) {
    seed = reducedSeedMgr_->instanceNew(seedName);
  }

  return seed.ptr();
}

/* Predicates */

struct HalfSliceIsInactive {
  explicit HalfSliceIsInactive(HalfSliceRank firstActive) :
    firstActive_(firstActive)
  {}

  bool operator()(const HalfTimeSlice::PtrConst & hs) const {
    return (hs->rank().value() - firstActive_.value() + (hs->direction() == HalfTimeSlice::BACKWARD)) % 2 != 0;
  }

private:
  HalfSliceRank firstActive_;
};

struct HalfSliceIsConverged {
  explicit HalfSliceIsConverged(HalfSliceRank firstActive) :
    firstActive_(firstActive)
  {}

  bool operator()(const HalfTimeSlice::PtrConst & hs) const {
    return hs->rank() < firstActive_;
  }

private:
  HalfSliceRank firstActive_;
};

/*template <typename T>
struct IsConverged {
  explicit IsConverged(HalfSliceRank firstActive) :
    firstActive_(firstActive)
  {}

  bool operator()(const typename std::map<HalfSliceRank, T>::value_type & v) const {
    return v.first < firstActive_;
  }

private:
  HalfSliceRank firstActive_;
};*/

void
LocalNetwork::convergedSlicesInc() {
  mapping_->convergedSlicesInc();
  halfTimeSlice_.erase(std::remove_if(halfTimeSlice_.begin(), halfTimeSlice_.end(),
                                      HalfSliceIsConverged(firstActiveSlice())),
                       halfTimeSlice_.end());

  leftSeedSync_.erase(leftSeedSync_.begin(), leftSeedSync_.lower_bound(firstActiveSlice()));
  jumpProjector_.erase(jumpProjector_.begin(), jumpProjector_.lower_bound(firstActiveSlice()));
  evenFullTimeSlice_.erase(evenFullTimeSlice_.begin(), evenFullTimeSlice_.lower_bound(firstActiveSlice()));
  oddFullTimeSlice_.erase(oddFullTimeSlice_.begin(), oddFullTimeSlice_.lower_bound(firstActiveSlice()));
  correctionSync_.erase(correctionSync_.begin(), correctionSync_.lower_bound(firstActiveSlice()));
  seedAssembler_.erase(seedAssembler_.begin(), seedAssembler_.lower_bound(firstActiveSlice()));

  evenCoarseTimeSlice_.erase(evenCoarseTimeSlice_.begin(), evenCoarseTimeSlice_.lower_bound(firstActiveSlice()));
  oddCoarseTimeSlice_.erase(oddCoarseTimeSlice_.begin(), oddCoarseTimeSlice_.lower_bound(firstActiveSlice()));
  fullCorrectionSync_.erase(fullCorrectionSync_.begin(), fullCorrectionSync_.lower_bound(firstActiveSlice()));
  
  mainSeed_.erase(mainSeed_.begin(), mainSeed_.lower_bound(firstActiveSlice()));
}

LocalNetwork::TaskList
LocalNetwork::halfTimeSlices() const {
  return TaskList(halfTimeSlice_.begin(), halfTimeSlice_.end());
}


LocalNetwork::TaskList
LocalNetwork::activeHalfTimeSlices() const {
  TaskList active;

  std::remove_copy_if(halfTimeSlice_.begin(), halfTimeSlice_.end(),
                      std::back_inserter(active),
                      HalfSliceIsInactive(firstActiveSlice()));

  return active;
}


LocalNetwork::TaskList
LocalNetwork::activeJumpProjectors() const {
  TaskList active;

  TaskMap::const_iterator it_end = jumpProjector_.end();
  for (TaskMap::const_iterator it = jumpProjector_.begin(); it != it_end; ++it) {
    if ((it->first.value() - mapping_->firstActiveSlice().value()) % 2 == 0) {
      active.push_back(it->second);
    }
  }

  return active;
}

LocalNetwork::TaskList
LocalNetwork::activeLeftSeedSyncs() const {
  TaskList active;

  TaskMap::const_iterator it_end = leftSeedSync_.end();
  for (TaskMap::const_iterator it = leftSeedSync_.begin(); it != it_end; ++it) {
    if ((it->first.value() - mapping_->firstActiveSlice().value()) % 2 == 0) {
      active.push_back(it->second);
    }
  }

  return active;
}

/*LocalNetwork::SRList
LocalNetwork::activeLeftSeedReaders() const {
  SRList active;

  SRMap::const_iterator it_end = leftSeedReader_.end();
  for (SRMap::const_iterator it = leftSeedReader_.begin(); it != it_end; ++it) {
    if ((it->first - mapping_->firstActiveSlice()) % 2 == 0) {
      active.push_back(it->second);
    }
  }

  return active;
}*/

/*
struct CorrectionTaskIsActive {
  explicit CorrectionTaskIsActive(HalfSliceRank firstActive) :
    firstActive_(firstActive)
   {}


  bool operator()(const LocalNetwork::TaskMap::value_type & p) const {
    return (p.first.value() - firstActive_.value()) % 2 == 0;
  }

private:
  HalfSliceRank firstActive_;
};*/


LocalNetwork::TaskList
LocalNetwork::activeFullTimeSlices() const {
  TaskList result;

  CpuRank firstActiveCpu = mapping_->hostCpu(mapping_->firstActiveSlice());
  CpuRank lastConvergedCpu = mapping_->hostCpu(mapping_->firstActiveSlice() - HalfSliceCount(1));

  if (localCpu() == firstActiveCpu || localCpu() == lastConvergedCpu) {
    ReducedSeed::Ptr firstActiveCorrection = reducedSeedMgr_->instance(toString(SeedId(SEED_CORRECTION, firstActiveSlice())));
    assert(firstActiveCorrection);
    
    UpdatedSeedAssembler::Ptr firstActiveAssembler = usaMgr_->instance(toString(firstActiveSlice()));
    assert(firstActiveAssembler);
    
    ZeroReducedCorrection::Ptr zeroCorrection = new ZeroReducedCorrection(firstActiveCorrection.ptr(), firstActiveAssembler.ptr());
    result.push_back(zeroCorrection);
  }

  const TaskMap & fullTimeSlice = ((firstActiveSlice().value()) % 2 == 0) ? evenFullTimeSlice_ : oddFullTimeSlice_;
  TaskMap::const_iterator it_end = fullTimeSlice.end();
  for (TaskMap::const_iterator it = fullTimeSlice.begin(); it != it_end; ++it) {
    result.push_back(it->second);
  }

  return result;
}

LocalNetwork::TaskList
LocalNetwork::activeCorrectionSyncs() const {
  TaskList result;

  TaskMap::const_iterator it_end = correctionSync_.end();
  for (TaskMap::const_iterator it = correctionSync_.begin(); it != it_end; ++it) {
    if ((it->first.value() - firstActiveSlice().value()) % 2 == 0) {
      result.push_back(it->second);
    }
  }

  return result;
}

LocalNetwork::TaskList
LocalNetwork::activeFullCorrectionSyncs() const {
  TaskList result;

  TaskMap::const_iterator it_end = fullCorrectionSync_.end();
  for (TaskMap::const_iterator it = fullCorrectionSync_.begin(); it != it_end; ++it) {
    if ((it->first.value() - firstActiveSlice().value()) % 2 == 0) {
      result.push_back(it->second);
    }
  }

  return result;
}

LocalNetwork::TaskList
LocalNetwork::activeCoarseTimeSlices() const {
  TaskList result;

  CpuRank firstActiveCpu = mapping_->hostCpu(mapping_->firstActiveSlice());
  CpuRank lastConvergedCpu = mapping_->hostCpu(mapping_->firstActiveSlice() - HalfSliceCount(1));

  if (localCpu() == firstActiveCpu || localCpu() == lastConvergedCpu) {
    Seed::Ptr firstActiveCorrection = seedMgr_->instance(toString(SeedId(SEED_CORRECTION, firstActiveSlice())));
    assert(firstActiveCorrection);
    
    UpdatedSeedAssembler::Ptr firstActiveAssembler = usaMgr_->instance(toString(firstActiveSlice()));
    assert(firstActiveAssembler);

    ZeroCorrection::Ptr zeroCorrection = new ZeroCorrection(firstActiveCorrection.ptr(), firstActiveAssembler.ptr());
    result.push_back(zeroCorrection);
  }

  const TaskMap & coarseTimeSlice = ((firstActiveSlice().value()) % 2 == 0) ? evenCoarseTimeSlice_ : oddCoarseTimeSlice_;
  TaskMap::const_iterator it_end = coarseTimeSlice.end();
  for (TaskMap::const_iterator it = coarseTimeSlice.begin(); it != it_end; ++it) {
    if ((it->first.value() - firstActiveSlice().value()) % 2 == 0) {
      result.push_back(it->second);
    }
  }

  return result;
}

LocalNetwork::TaskList
LocalNetwork::activeSeedAssemblers() const {
  TaskList result;

  TaskMap::const_iterator it_end = seedAssembler_.end();
  for (TaskMap::const_iterator it = seedAssembler_.begin(); it != it_end; ++it) {
    if ((it->first.value() - firstActiveSlice().value()) % 2 == 0) {
      result.push_back(it->second);
    }
  }

  return result;
}


LocalNetwork::SeedMap
LocalNetwork::mainSeeds() const {
  return SeedMap(mainSeed_);
}

LocalNetwork::MainSeedMap
LocalNetwork::activeMainSeeds() const {
  MainSeedMap msp;

  SeedMap::const_iterator it_end = mainSeed_.end();
  for (SeedMap::const_iterator it = mainSeed_.begin(); it != it_end; ++it) {
    if ((it->first.value() - mapping_->firstActiveSlice().value()) % 2 == 0) {
      msp.insert(msp.end(), std::make_pair(SliceRank(it->first.value() / 2), it->second));
    }
  }

  return msp;
}

} /* end namespace Hts */ } /* end namespace Pita */
