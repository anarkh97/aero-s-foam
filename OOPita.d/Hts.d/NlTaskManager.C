#include "NlTaskManager.h"

#include "BasicHalfTimeSliceImpl.h"
#include "../JumpBuilder.h"
#include "../ConcurrentBasisManager.h"
#include "../NearSymmetricSolver.h"
#include "../DynamStateReductor.h"
#include "../DynamStateReconstructor.h"
#include "NlProjectionNetwork.h"
#include "PropagationDataSharing.h"

#include "../InitialSeedTask.h"

namespace Pita { namespace Hts {

NlTaskManager::NlTaskManager(SliceMapping * mapping, RemoteState::MpiManager * commMgr,
                             NlPropagatorManager * propagatorMgr, SeedInitializer * seedInitializer,
                             PostProcessing::Manager * postProcessingMgr, SeedErrorEvaluator::Manager * jumpEvaluatorMgr,
                             double projectorTolerance, IterationRank lastIteration) :
  TaskManager(IterationRank(0)),
  mapping_(mapping),
  commMgr_(commMgr),
  propagatorMgr_(propagatorMgr),
  seedInitializer_(seedInitializer),
  postProcessingMgr_(postProcessingMgr),
  jumpEvaluatorMgr_(jumpEvaluatorMgr),
  projectorTolerance_(projectorTolerance),
  phase_(NULL),
  continuation_(&NlTaskManager::noop),
  localNetwork_(NULL),
  sharing_(NULL),
  projectionNetwork_(NULL),
  lastIteration_(lastIteration)
{
  initialize();
  scheduleInitialSeed();
}

void
NlTaskManager::iterationInc() {
  setIteration(iteration().next());
  scheduleDataSharing();
}

NlTaskManager::Phase *
NlTaskManager::phase() {
  return phase_.ptr();
}

void
NlTaskManager::phaseInc() {
  (this->*continuation_)();
}

void
NlTaskManager::initialize() {
  // Projection (global data sharing)
  sharing_ = new PropagationDataSharing(mapping_.ptr(), commMgr_->localCpu(), commMgr_->communicator(), commMgr_->vectorSize());
  NlDynamOps::Ptr dynOps = propagatorMgr_->dynamOpsNew();
  projectionNetwork_ = new NlProjectionNetwork(sharing_.ptr(), dynOps.ptr(), commMgr_->vectorSize(), projectorTolerance_);

  // Parallel time-integration
  propagatorMgr_->postProcessingManagerIs(postProcessingMgr_.ptr());
  propagatorMgr_->concurrentBasisManagerIs(projectionNetwork_->concurrentMgr());
  propagatorMgr_->propagatedBasisManagerIs(projectionNetwork_->endBasisMgr());
  propagatorMgr_->collectorIs(sharing_->collector());
  BasicHalfTimeSliceImpl::Manager::Ptr htsMgr = BasicHalfTimeSliceImpl::Manager::New(propagatorMgr_.ptr());

  // Recurring tasks 
  localNetwork_ = new NlLocalNetwork(mapping_.ptr(),
                                     commMgr_.ptr(),
                                     htsMgr.ptr(),
                                     projectionNetwork_->corrRedMgr(),
                                     projectionNetwork_->corrReconMgr(),
                                     projectionNetwork_->condensMgr(),
                                     projectionNetwork_->projBuildMgr(),
                                     jumpEvaluatorMgr_.ptr());
  localNetwork_->statusIs(LocalNetwork::ACTIVE);
}


void
NlTaskManager::scheduleNothing() {
  setPhase(NULL);
  setContinuation(&NlTaskManager::noop);
}

void
NlTaskManager::scheduleInitialSeed() {
  TaskList initialSeedInformation;

  typedef LocalNetwork::MainSeedMap MainSeedMap;
  MainSeedMap primalSeeds = localNetwork_->activeSeeds();
  for (MainSeedMap::iterator it = primalSeeds.begin(); it != primalSeeds.end(); ++it) {
    Seed::Ptr targetSeed = it->second;
    Seed::Status initialSeedStatus = (it->first == SliceRank(0)) ? Seed::CONVERGED : Seed::ACTIVE;
    NamedTask::Ptr task = new InitialSeedTask(targetSeed.ptr(), seedInitializer_.ptr(), it->first, initialSeedStatus);
    initialSeedInformation.push_back(task);
  }
  
  setPhase(phaseNew("Initial Seed Information", initialSeedInformation));
  setContinuation(&NlTaskManager::fillFirstProjectionBasis); 
}

void
NlTaskManager::fillFirstProjectionBasis() {
  int initSeedCount = (mapping_->activeSlices().value() / 2) + 1;
  for (int i = 0; i < initSeedCount; ++i) {
    projectionNetwork_->sharedProjectionBasis()->lastStateIs(seedInitializer_->initialSeed(SliceRank(i)));
  }

  scheduleProjectionBasisCondensation();
}

void
NlTaskManager::scheduleProjectionBasisCondensation() {
  LocalNetwork::TaskList oldStyleList = localNetwork_->activeCondensations();
  setPhase(phaseNew("Projection Basis Condensation", TaskList(oldStyleList.begin(), oldStyleList.end())));
  setContinuation(&NlTaskManager::scheduleFinePropagation);
}

void
NlTaskManager::scheduleFinePropagation() {
  LocalNetwork::TaskList oldStyleList = localNetwork_->activeFinePropagators();
  setPhase(phaseNew("Parallel Fine Propagation", TaskList(oldStyleList.begin(), oldStyleList.end())));
  setContinuation(&NlTaskManager::scheduleNothing);
}

void
NlTaskManager::scheduleDataSharing() {
  TaskList dataSharing;
  if (iteration() < lastIteration_) {
    dataSharing.push_back(sharing_);
  }
  setPhase(phaseNew("Data Sharing", dataSharing));
  setContinuation(&NlTaskManager::checkConvergence);
}

void
NlTaskManager::checkConvergence() {
  localNetwork_->convergedSlicesInc();
  schedulePropagatedSeedSynchronization();
}

void
NlTaskManager::schedulePropagatedSeedSynchronization() {
  LocalNetwork::TaskList oldStyleList = localNetwork_->activePropagatedSeedSyncs();
  setPhase(phaseNew("Propagated Seed Synchronization", TaskList(oldStyleList.begin(), oldStyleList.end())));
  setContinuation(&NlTaskManager::scheduleJumpProjection);
}

void
NlTaskManager::scheduleJumpProjection() {
  LocalNetwork::TaskList oldStyleList = localNetwork_->activeJumpBuilders();
  setPhase(phaseNew("Jump Evaluation", TaskList(oldStyleList.begin(), oldStyleList.end())));
  setContinuation(&NlTaskManager::scheduleProjectionBuilding);
}

void
NlTaskManager::scheduleProjectionBuilding() {
  LocalNetwork::TaskList oldStyleList = localNetwork_->activeProjectionBuilders();
  setPhase(phaseNew("Projection Building", TaskList(oldStyleList.begin(), oldStyleList.end())));
  setContinuation(&NlTaskManager::scheduleCorrectionPropagation);
}

void
NlTaskManager::scheduleCorrectionPropagation() {
  LocalNetwork::TaskList oldStyleList = localNetwork_->activeCorrectionPropagators();
  setPhase(phaseNew("Correction Propagation", TaskList(oldStyleList.begin(), oldStyleList.end())));
  setContinuation(&NlTaskManager::scheduleSeedUpdate);
}

void
NlTaskManager::scheduleSeedUpdate() {
  LocalNetwork::TaskList oldStyleList = localNetwork_->activeSeedUpdaters();
  setPhase(phaseNew("Seed Update", TaskList(oldStyleList.begin(), oldStyleList.end())));
  setContinuation(&NlTaskManager::enrichProjectionBasis);
}

void
NlTaskManager::enrichProjectionBasis() {
  projectionNetwork_->sharedProjectionBasis()->lastStateBasisIs(sharing_->consolidatedBasis());
  if (iteration() == lastIteration_) { 
    projectionNetwork_->sharedProjectionBasis()->stateBasisDel();
  }
  scheduleProjectionBasisCondensation();
}

} /* end namespace Hts */ } /* end namespace Pita */
