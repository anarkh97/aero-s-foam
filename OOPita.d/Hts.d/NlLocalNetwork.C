#include "NlLocalNetwork.h"

#include "../ReducedSeedWriterTask.h"

namespace Pita { namespace Hts {

NlLocalNetwork::NlLocalNetwork(SliceMapping * mapping,
                               RemoteState::MpiManager * commMgr,
                               HalfTimeSlice::Manager * htsMgr,
                               CorrectionReductor::Manager * corrRedMgr,
                               CorrectionReconstructor::Manager * corrReconMgr,
                               BasisCondensationManager * condensMgr,
                               ProjectionBuildingFactory * projBuildMgr,
                               SeedErrorEvaluator::Manager * jumpEvalMgr) :
  LocalNetwork(mapping, commMgr),
  htsMgr_(htsMgr),
  jumpMgr_(JumpBuilder::ManagerImpl::New()),
  seedUpMgr_(SeedUpdater::ManagerImpl::New()),
  corrRedMgr_(corrRedMgr),
  corrReconMgr_(corrReconMgr),
  condensMgr_(condensMgr),
  projBuildMgr_(projBuildMgr),
  jumpEvalMgr_(jumpEvalMgr),
  noCorrectionMgr_(NoCorrectionManager::New())
{}

void
NlLocalNetwork::init() {
  for (SliceMapping::SliceIterator it = hostedSlice(localCpu()); it; ++it) {
    HalfSliceRank sliceRank = *it;
    
    HalfSliceRank beginSeedRank = sliceRank;
    HalfSliceRank endSeedRank = sliceRank.next();

    addPropagatedSeedSend(endSeedRank);
    addCorrectionReconstructor(endSeedRank);
    addReducedCorrectionRecv(endSeedRank);
    addSeedUpdater(endSeedRank);
    addCorrectionSend(endSeedRank);
    
    addNoCorrection(endSeedRank);
    
    addMainSeed(endSeedRank);
    addBackwardCondensation(sliceRank);
    addBackwardPropagation(sliceRank);

    addPropagatedSeedRecv(beginSeedRank);
    addJumpBuilder(beginSeedRank);
    addCorrectionRecv(beginSeedRank);
    addCorrectionReductor(beginSeedRank);
    addProjectionBuilding(beginSeedRank);
    addReducedCorrectionSend(beginSeedRank + FullSliceCount(1));
    addSeedUpdater(beginSeedRank);
    
    addMainSeed(beginSeedRank);
    addForwardCondensation(sliceRank);
    addForwardPropagation(sliceRank);

    addNoCorrection(beginSeedRank);
  }
}

void
NlLocalNetwork::addForwardPropagation(HalfSliceRank sliceRank) {
  NamedTask::Ptr task = forwardHalfSliceNew(sliceRank);

  int iterParity = parity(sliceRank);
  finePropagators_[iterParity][sliceRank] = task;
}

void
NlLocalNetwork::addBackwardPropagation(HalfSliceRank sliceRank) {
  if (sliceRank == HalfSliceRank(0)) return; // Domain guard

  NamedTask::Ptr task = backwardHalfSliceNew(sliceRank);

  int iterParity = parity(sliceRank.next());
  finePropagators_[iterParity][sliceRank] = task;
}

void
NlLocalNetwork::addForwardCondensation(HalfSliceRank sliceRank) {
  NamedTask::Ptr task = condensMgr_->instanceNew(HalfSliceId(sliceRank, FORWARD));

  int iterParity = parity(sliceRank);
  condensations_[iterParity][sliceRank] = task;
}

void
NlLocalNetwork::addBackwardCondensation(HalfSliceRank sliceRank) {
  if (sliceRank == HalfSliceRank(0)) return; // Domain guard
  
  NamedTask::Ptr task = condensMgr_->instanceNew(HalfSliceId(sliceRank, BACKWARD));

  int iterParity = parity(sliceRank.previous());
  condensations_[iterParity][sliceRank] = task;
}

void
NlLocalNetwork::addMainSeed(HalfSliceRank seedRank) {
  int iterParity = parity(seedRank);
  SliceRank fullSeedRank(seedRank.value() / 2); 
  seeds_[iterParity][fullSeedRank] = fullSeedGet(SeedId(MAIN_SEED, seedRank));
}

void
NlLocalNetwork::addPropagatedSeedSend(HalfSliceRank seedRank) {
  if (seedRank >= firstInactiveSlice()) return; // Domain guard

  NamedTask::Ptr task = propagatedSeedSendNew(seedRank);
  if (!task) return; // Comm guard

  int iterParity = parity(seedRank);
  propagatedSeedSyncs_[iterParity][seedRank] = task;
}

void
NlLocalNetwork::addPropagatedSeedRecv(HalfSliceRank seedRank) {
  if (seedRank == HalfSliceRank(0)) return; // Domain guard

  NamedTask::Ptr task = propagatedSeedRecvNew(seedRank);
  if (!task) return; // Comm guard

  int iterParity = parity(seedRank);
  propagatedSeedSyncs_[iterParity][seedRank] = task;
}

void
NlLocalNetwork::addJumpBuilder(HalfSliceRank seedRank) {
  if (seedRank == HalfSliceRank(0)) return; // Domain guard

  NamedTask::Ptr task = jumpBuilderNew(seedRank);

  int iterParity = parity(seedRank);
  jumpBuilders_[iterParity][seedRank] = task;
  
  if (jumpEvalMgr()) { 
    SeedErrorEvaluator::Ptr jumpEvaluator = jumpEvalMgr()->instanceNew(fullSeedGet(SeedId(SEED_JUMP, seedRank)));
    jumpEvaluator->referenceSeedIs(fullSeedGet(SeedId(LEFT_SEED, seedRank)));
  }
}

void
NlLocalNetwork::addProjectionBuilding(HalfSliceRank seedRank) {
  if (seedRank == HalfSliceRank(0)) return; // Domain guard
  if (seedRank + HalfSliceCount(2) > firstInactiveSlice()) return; // Domain guard

  NamedTask::Ptr task = projBuildMgr_->instanceNew(seedRank);

  int iterParity = parity(seedRank);
  projectionBuilders_[iterParity][seedRank] = task;
}

void
NlLocalNetwork::addCorrectionReductor(HalfSliceRank seedRank) {
  if (seedRank == HalfSliceRank(0)) return; // Domain guard
  if (seedRank + HalfSliceCount(2) > firstInactiveSlice()) return; // Domain guard

  NamedTask::Ptr task = correctionReductorNew(seedRank);

  int iterParity = parity(seedRank);
  correctionPropagators_[iterParity][ActivationRange(seedRank, HalfSliceCount(1))] = task; // TODO Ordering relative to reconstructor
}

void
NlLocalNetwork::addCorrectionReconstructor(HalfSliceRank seedRank) {
  if (seedRank > firstInactiveSlice()) return; // Domain guard
  if (seedRank <= HalfSliceRank(2)) return;

  NamedTask::Ptr task = correctionReconstructorNew(seedRank);

  int iterParity = parity(seedRank);
  correctionPropagators_[iterParity][ActivationRange(seedRank.previous(), HalfSliceCount(-1))] = task; // TODO Ordering relative to reductor 
}

void
NlLocalNetwork::addCorrectionSend(HalfSliceRank seedRank) {
  if (seedRank > firstInactiveSlice()) return; // Domain guard
  if (seedRank <= HalfSliceRank(2)) return;

  NamedTask::Ptr task = correctionSendNew(seedRank);
  if (!task) return; // Comm guard

  int iterParity = parity(seedRank);
  correctionPropagators_[iterParity][ActivationRange(seedRank, HalfSliceCount(-1))] = task; // TODO Ordering relative to reconstructor/reductor
}

void
NlLocalNetwork::addCorrectionRecv(HalfSliceRank seedRank) {
  if (seedRank > firstInactiveSlice()) return; // Domain guard
  if (seedRank <= HalfSliceRank(2)) return; // Domain guard

  NamedTask::Ptr task = correctionRecvNew(seedRank);
  if (!task) return; // Comm guard

  int iterParity = parity(seedRank);
  correctionPropagators_[iterParity][ActivationRange(seedRank.previous(), HalfSliceCount(1))] = task; // TODO Ordering relative to reconstructor/reductor
}

void
NlLocalNetwork::addReducedCorrectionSend(HalfSliceRank seedRank) {
  if (seedRank > firstInactiveSlice()) return; // Domain guard
  if (seedRank <= HalfSliceRank(2)) return; // Domain guard

  NamedTask::Ptr task = reducedCorrectionSendNew(seedRank);
  if (!task) return; // Comm guard

  int iterParity = parity(seedRank);
  correctionPropagators_[iterParity][ActivationRange(seedRank.previous(), HalfSliceCount(-1))] = task; // TODO Ordering relative to reconstructor/reductor
}

void
NlLocalNetwork::addReducedCorrectionRecv(HalfSliceRank seedRank) {
  if (seedRank > firstInactiveSlice()) return; // Domain guard
  if (seedRank <= HalfSliceRank(2)) return; // Domain guard

  NamedTask::Ptr task = reducedCorrectionRecvNew(seedRank);
  if (!task) return; // Comm guard

  int iterParity = parity(seedRank);
  correctionPropagators_[iterParity][ActivationRange(seedRank - FullSliceCount(1), HalfSliceCount(1))] = task; // TODO Ordering relative to reconstructor/reductor
}

void
NlLocalNetwork::addSeedUpdater(HalfSliceRank seedRank) {
  if (seedRank == HalfSliceRank(0)) return; // Domain guard

  NamedTask::Ptr task = seedUpdater(seedRank);
  if (!task) {
    task = seedUpdaterNew(seedRank);
  }

  int iterParity = parity(seedRank);
  seedUpdaters_[iterParity][seedRank] = task;
}

void
NlLocalNetwork::addNoCorrection(HalfSliceRank seedRank) {
  noCorrectionMgr_->notifierIs(fullSeedGet(SeedId(SEED_CORRECTION, seedRank)), fullSeedGet(SeedId(LEFT_SEED, seedRank)));
}

void
NlLocalNetwork::statusIs(Status s) {
  if (status() == s) return;

  if (s == ACTIVE) {
    init();
    return;
  } 

  throw Fwk::RangeException();  // Not implemented
}

void
NlLocalNetwork::convergedSlicesInc() {
  LocalNetwork::convergedSlicesInc();

  for (int parity = 0; parity < 2; ++parity) {
    HalfSliceRank start = firstActiveSlice();
    finePropagators_[parity].erase(finePropagators_[parity].begin(), finePropagators_[parity].lower_bound(start));
    propagatedSeedSyncs_[parity].erase(propagatedSeedSyncs_[parity].begin(), propagatedSeedSyncs_[parity].lower_bound(start));
    jumpBuilders_[parity].erase(jumpBuilders_[parity].begin(), jumpBuilders_[parity].lower_bound(start));
    correctionPropagators_[parity].erase(correctionPropagators_[parity].begin(), correctionPropagators_[parity].lower_bound(start));
    seedUpdaters_[parity].erase(seedUpdaters_[parity].begin(), seedUpdaters_[parity].lower_bound(start));

    condensations_[parity].erase(condensations_[parity].begin(), condensations_[parity].lower_bound(start));
    projectionBuilders_[parity].erase(projectionBuilders_[parity].begin(), projectionBuilders_[parity].lower_bound(start));

    SliceRank seedStart((start.value() + parity) / 2);
    seeds_[parity].erase(seeds_[parity].begin(), seeds_[parity].lower_bound(seedStart));
  }
}

NamedTask::Ptr
NlLocalNetwork::forwardHalfSliceNew(HalfSliceRank sliceRank) {
  HalfTimeSlice::Ptr result = htsMgr()->instanceNew(HalfSliceId(sliceRank, FORWARD));

  result->seedIs(fullSeedGet(SeedId(MAIN_SEED, sliceRank)));
  result->propagatedSeedIs(fullSeedGet(SeedId(LEFT_SEED, sliceRank.next())));

  return result;
}

NamedTask::Ptr
NlLocalNetwork::backwardHalfSliceNew(HalfSliceRank sliceRank) {
  HalfTimeSlice::Ptr result = htsMgr()->instanceNew(HalfSliceId(sliceRank, BACKWARD));
  
  result->seedIs(fullSeedGet(SeedId(MAIN_SEED, sliceRank.next())));
  result->propagatedSeedIs(fullSeedGet(SeedId(RIGHT_SEED, sliceRank)));

  return result;
}

NamedTask::Ptr
NlLocalNetwork::propagatedSeedSendNew(HalfSliceRank seedRank) {
  assert(hostCpu(seedRank.previous()) == localCpu());
  CpuRank targetCpu(hostCpu(seedRank));
  
  RemoteStateTask::Ptr result = NULL;
  if (targetCpu != localCpu()) {
    RemoteState::SeedReader::Ptr reader = commMgr()->readerNew(fullSeedGet(SeedId(LEFT_SEED, seedRank)), targetCpu);
    String taskName = String("Send Propagated Seed ") + toString(seedRank);
    result = RemoteStateTask::New(taskName, reader.ptr());
  }
  return result;
}

NamedTask::Ptr
NlLocalNetwork::propagatedSeedRecvNew(HalfSliceRank seedRank) {
  assert(hostCpu(seedRank) == localCpu());
  CpuRank previousCpu = hostCpu(seedRank.previous());

  RemoteStateTask::Ptr result = NULL;
  if (previousCpu != localCpu()) {
    RemoteState::SeedWriter::Ptr writer = commMgr()->writerNew(fullSeedGet(SeedId(LEFT_SEED, seedRank)), previousCpu);
    String taskName = String("Receive Propagated Seed ") + toString(seedRank);
    result = RemoteStateTask::New(taskName, writer.ptr());
  }
  return result;
}

NamedTask::Ptr
NlLocalNetwork::jumpBuilderNew(HalfSliceRank seedRank) {
  JumpBuilder::Ptr result = jumpMgr()->instanceNew(toString(seedRank));

  result->predictedSeedIs(fullSeedGet(SeedId(RIGHT_SEED, seedRank)));
  result->actualSeedIs(fullSeedGet(SeedId(LEFT_SEED, seedRank)));
  result->seedJumpIs(fullSeedGet(SeedId(SEED_JUMP, seedRank)));

  return result;
}

NamedTask::Ptr
NlLocalNetwork::correctionReductorNew(HalfSliceRank seedRank) {
  CorrectionReductor::Ptr result = corrRedMgr()->instanceNew(seedRank);

  result->jumpIs(fullSeedGet(SeedId(SEED_JUMP, seedRank)));
  result->correctionIs(fullSeedGet(SeedId(SEED_CORRECTION, seedRank)));
  result->nextCorrectionIs(reducedSeedGet(SeedId(SEED_CORRECTION, seedRank + HalfSliceCount(2))));

  return result;
}

NamedTask::Ptr
NlLocalNetwork::correctionReconstructorNew(HalfSliceRank seedRank) {
  CorrectionReconstructor::Ptr result = corrReconMgr()->instanceNew(seedRank);

  result->correctionIs(fullSeedGet(SeedId(SEED_CORRECTION, seedRank)));
  result->correctionComponentsIs(reducedSeedGet(SeedId(SEED_CORRECTION, seedRank)));

  return result;
}

NamedTask::Ptr
NlLocalNetwork::correctionSendNew(HalfSliceRank seedRank) {
  CpuRank nextHeadCpu = hostCpu(seedRank);

  RemoteStateTask::Ptr result = NULL;
  if (nextHeadCpu != localCpu()) {
    RemoteState::SeedReader::Ptr reader = commMgr()->readerNew(fullSeedGet(SeedId(SEED_CORRECTION, seedRank)), nextHeadCpu);
    String taskName = String("Send Correction ") + toString(seedRank);
    result = RemoteStateTask::New(taskName, reader.ptr());
  }

  return result;
}

NamedTask::Ptr
NlLocalNetwork::correctionRecvNew(HalfSliceRank seedRank) {
  CpuRank previousTailCpu = hostCpu(seedRank.previous());

  RemoteStateTask::Ptr result = NULL;
  if (previousTailCpu != localCpu()) {
    RemoteState::SeedWriter::Ptr writer = commMgr()->writerNew(fullSeedGet(SeedId(SEED_CORRECTION, seedRank)), previousTailCpu);
    String taskName = String("Recv Correction ") + toString(seedRank);
    result = RemoteStateTask::New(taskName, writer.ptr());
  }

  return result;
}

NamedTask::Ptr
NlLocalNetwork::reducedCorrectionSendNew(HalfSliceRank seedRank) {
  CpuRank tailCpu = hostCpu(seedRank.previous());

  RemoteStateTask::Ptr result = NULL;
  if (tailCpu != localCpu()) {
    RemoteState::ReducedSeedReader::Ptr reader = commMgr()->readerNew(reducedSeedGet(SeedId(SEED_CORRECTION, seedRank)), tailCpu);
    String taskName = String("Send Reduced Correction ") + toString(seedRank);
    result = RemoteStateTask::New(taskName, reader.ptr());
  }

  return result;
}

NamedTask::Ptr
NlLocalNetwork::reducedCorrectionRecvNew(HalfSliceRank seedRank) {
  CpuRank headCpu = hostCpu(seedRank - FullSliceCount(1));

  RemoteStateTask::Ptr result = NULL;
  if (headCpu != localCpu()) {
    RemoteState::ReducedSeedWriter::Ptr writer = commMgr()->writerNew(reducedSeedGet(SeedId(SEED_CORRECTION, seedRank)), headCpu);
    RemoteState::MpiReducedSeedWriter::Ptr mpiWriter = ptr_cast<RemoteState::MpiReducedSeedWriter>(writer); // Safe cast
    String taskName = String("Recv Reduced Correction ") + toString(seedRank);
    CorrectionReconstructor::Ptr reconstructor = corrReconMgr()->instance(seedRank);
    assert(reconstructor);
    result = ReducedSeedWriterTask<CorrectionReconstructor>::New(taskName, mpiWriter.ptr(), reconstructor.ptr());
  }

  return result;
}

NamedTask::Ptr
NlLocalNetwork::seedUpdater(HalfSliceRank seedRank) {
  return seedUpMgr()->instance(toString(seedRank));
}

NamedTask::Ptr
NlLocalNetwork::seedUpdaterNew(HalfSliceRank seedRank) {
  SeedUpdater::Ptr result = seedUpMgr()->instanceNew(toString(seedRank));

  result->propagatedSeedIs(fullSeedGet(SeedId(LEFT_SEED, seedRank)));
  result->correctionIs(fullSeedGet(SeedId(SEED_CORRECTION, seedRank)));
  result->updatedSeedIs(fullSeedGet(SeedId(MAIN_SEED, seedRank)));

  return result;
}

} /* end namespace Hts */ } /* end namespace Pita */
