#include "LocalNetwork.h"

#include "RemoteStateTask.h"

#include <algorithm>
#include <iterator>

namespace Pita {
namespace Hts {

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

LocalNetwork::LocalNetwork(FullSliceCount totalFullSlices,
                           CpuCount availableCpus,
                           HalfSliceCount maxWorkload,
                           CpuRank localCpu,
                           HalfTimeSlice::Manager * sliceMgr,
                           JumpProjector::Manager * jumpProjMgr,
                           ReducedFullTimeSlice::Manager * ftsMgr,
                           RemoteState::Manager * commMgr) :
  totalFullSlices_(totalFullSlices),
  availableCpus_(availableCpus),
  maxWorkload_(maxWorkload),
  localCpu_(localCpu),
  taskManager_(TaskManager::New(2 * totalFullSlices.value(), availableCpus.value(), maxWorkload.value())),
  sliceMgr_(sliceMgr),
  jumpProjMgr_(jumpProjMgr),
  ftsMgr_(ftsMgr),
  commMgr_(commMgr),
  seedMgr_(Seed::Manager::New()),
  reducedSeedMgr_(ReducedSeed::Manager::New())
{
  init();
}

void
LocalNetwork::init() {
  for (TaskManager::TaskIterator it = taskManager_->tasks(localCpu().value()); it; ++it) {
    HalfSliceRank sliceRank(*it);

    // HalfTimeSlices
    HalfTimeSlice::Ptr forwardSlice = sliceMgr_->instanceNew(HalfSliceId(sliceRank, HalfTimeSlice::FORWARD));
    forwardSlice->seedIs(getSeed(SeedId(MAIN_SEED, sliceRank)));
    Seed::Ptr nextLeftPropagatedSeed = getSeed(SeedId(LEFT_SEED, sliceRank + HalfSliceCount(1)));
    forwardSlice->propagatedSeedIs(nextLeftPropagatedSeed.ptr());
    halfTimeSlice_.push_back(forwardSlice);

    HalfTimeSlice::Ptr backwardSlice = sliceMgr_->instanceNew(HalfSliceId(sliceRank, HalfTimeSlice::BACKWARD));
    backwardSlice->seedIs(getSeed(SeedId(MAIN_SEED, sliceRank + HalfSliceCount(1))));
    backwardSlice->propagatedSeedIs(getSeed(SeedId(RIGHT_SEED, sliceRank)));
    halfTimeSlice_.push_back(backwardSlice);

    // JumpProjectors
    JumpProjector::Ptr jumpProjector = jumpProjMgr_->instanceNew(toString(sliceRank));
    jumpProjector->predictedSeedIs(getSeed(SeedId(RIGHT_SEED, sliceRank)));
    Seed::Ptr leftPropagatedSeed = getSeed(SeedId(LEFT_SEED, sliceRank));
    jumpProjector->actualSeedIs(leftPropagatedSeed.ptr());
    jumpProjector->seedJumpIs(getSeed(SeedId(SEED_JUMP, sliceRank)));
    jumpProjector->reducedSeedJumpIs(getReducedSeed(SeedId(SEED_JUMP, sliceRank)));
    jumpProjector_.insert(jumpProjector_.end(), std::make_pair(sliceRank, jumpProjector));

    // FullTimeSlices
    if (sliceRank  + HalfSliceCount(1) < firstInactiveSlice()) {
      ReducedFullTimeSlice::Ptr fullSlice = ftsMgr_->instanceNew(sliceRank);
      fullSlice->jumpIs(getReducedSeed(SeedId(SEED_JUMP, sliceRank)));
      ReducedSeed::Ptr correction = getReducedSeed(SeedId(SEED_CORRECTION, sliceRank));
      fullSlice->correctionIs(correction.ptr());
      ReducedSeed::Ptr nextCorrection = getReducedSeed(SeedId(SEED_CORRECTION, sliceRank + HalfSliceCount(2)));
      fullSlice->nextCorrectionIs(nextCorrection.ptr());
      //log() << "Setup FTS head/tail = " << fullSlice->headHalfSlice() << "/" << fullSlice->tailHalfSlice() << "\n";

      HalfSliceRank previousSliceRank = sliceRank - HalfSliceCount(2);
      CpuRank previousCpu(taskManager_->worker(previousSliceRank.value()));
      if (previousCpu != CpuRank(-1) && previousCpu != localCpu()) {
        RemoteState::ReducedSeedWriter::Ptr writer = commMgr_->writerNew(correction.ptr(), previousCpu);
        if (writer) {
          String taskName = String("Receive Correction ") + toString(previousSliceRank);
          RemoteStateTask::Ptr task = RemoteStateTask::New(taskName, writer.ptr());
          fullTimeSlice_.insert(fullTimeSlice_.end(), std::make_pair(previousSliceRank, task));
        }
      }

      fullTimeSlice_.insert(fullTimeSlice_.end(), std::make_pair(sliceRank, fullSlice));
      
      HalfSliceRank nextSliceRank = sliceRank + HalfSliceCount(2);
      CpuRank nextCpu(taskManager_->worker(nextSliceRank.value()));
      if (nextCpu != CpuRank(-1) && nextCpu != localCpu()) {
        RemoteState::ReducedSeedReader::Ptr reader = commMgr_->readerNew(nextCorrection.ptr(), nextCpu);
        if (reader) {
          String taskName = String("Send Correction ") + toString(nextSliceRank);
          RemoteStateTask::Ptr task = RemoteStateTask::New(taskName, reader.ptr());
          fullTimeSlice_.insert(fullTimeSlice_.end(), std::make_pair(nextSliceRank, task));
        }
      }
    }

    CpuRank writerCpu(taskManager_->worker(sliceRank.value() - 1));
    if (writerCpu != CpuRank(-1) && writerCpu != localCpu()) {
      //log() << "Setup writer from Cpu " << writerCpu << " for " << leftPropagatedSeed->name()
      //  << " and " << "slice " << sliceRank - HalfSliceCount(1) << "\n";
      RemoteState::SeedWriter::Ptr writer = commMgr_->writerNew(leftPropagatedSeed.ptr(), writerCpu);
      if (writer) {
        leftSeedWriter_.insert(leftSeedWriter_.end(),
                               std::make_pair(sliceRank - HalfSliceCount(1), writer));
      }
    }

    CpuRank readerCpu(taskManager_->worker(sliceRank.value() + 1));
    if (readerCpu != CpuRank(-1) && readerCpu != localCpu()) {
      //log() << "Setup reader to Cpu " << readerCpu << " for " << nextLeftPropagatedSeed->name() <<
      //  " and " << "slice " << sliceRank << "\n";
      RemoteState::SeedReader::Ptr reader = commMgr_->readerNew(nextLeftPropagatedSeed.ptr(), readerCpu);
      if (reader) {
        leftSeedReader_.insert(leftSeedReader_.end(),
                               std::make_pair(sliceRank, reader));
      }
    }

    // Update
  
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

void
LocalNetwork::convergedSlicesInc() {
  taskManager_->completedTasksInc(1);
  halfTimeSlice_.erase(std::remove_if(halfTimeSlice_.begin(), halfTimeSlice_.end(),
                                      HalfSliceIsConverged(firstActiveSlice())),
                       halfTimeSlice_.end());
  // TODO other tasks
}

LocalNetwork::HTSList
LocalNetwork::activeHalfTimeSlices() const {
  HTSList active;

  std::remove_copy_if(halfTimeSlice_.begin(), halfTimeSlice_.end(),
                      std::back_inserter(active),
                      HalfSliceIsInactive(firstActiveSlice()));

  return active;
}


LocalNetwork::JPList
LocalNetwork::activeJumpProjectors() const {
  JPList active;

  JPMap::const_iterator it_end = jumpProjector_.end();
  for (JPMap::const_iterator it = jumpProjector_.begin(); it != it_end; ++it) {
    if ((it->first.value() - taskManager_->firstCurrentTask()) % 2 == 0) {
      active.push_back(it->second);
    }
  }

  return active;
}

LocalNetwork::SWList
LocalNetwork::activeLeftSeedWriters() const {
  SWList active;

  SWMap::const_iterator it_end = leftSeedWriter_.end();
  for (SWMap::const_iterator it = leftSeedWriter_.begin(); it != it_end; ++it) {
    if ((it->first.value() - taskManager_->firstCurrentTask()) % 2 == 0) {
      active.push_back(it->second);
    }
  }

  return active;
}

LocalNetwork::SRList
LocalNetwork::activeLeftSeedReaders() const {
  SRList active;

  SRMap::const_iterator it_end = leftSeedReader_.end();
  for (SRMap::const_iterator it = leftSeedReader_.begin(); it != it_end; ++it) {
    if ((it->first.value() - taskManager_->firstCurrentTask()) % 2 == 0) {
      active.push_back(it->second);
    }
  }

  return active;
}
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
 
  TaskMap::const_iterator it_end = fullTimeSlice_.end();
  for (TaskMap::const_iterator it = fullTimeSlice_.begin(); it != it_end; ++it) {
    if ((it->first.value() - firstActiveSlice().value()) % 2 == 0) {
      result.push_back(it->second);
    }
  }

  return result;
}

LocalNetwork::MainSeedMap
LocalNetwork::activeMainSeeds() const {
  MainSeedMap msp;

  SeedMap::const_iterator it_end = mainSeed_.end();
  for (SeedMap::const_iterator it = mainSeed_.begin(); it != it_end; ++it) {
    if ((it->first.value() - taskManager_->firstCurrentTask()) % 2 == 0) {
      msp.insert(msp.end(), std::make_pair(SliceRank(it->first.value() / 2), it->second));
    }
  }

  return msp;
}

} /* end namespace Hts */ } /* end namespace Pita */
