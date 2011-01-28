#include "NlTaskManager.h"

#include "../JumpBuilder.h"
#include "../ConcurrentBasisManager.h"
#include "../NearSymmetricSolver.h"
#include "../DynamStateReductor.h"
#include "../DynamStateReconstructor.h"

#include "../InitialSeedTask.h"

namespace Pita { namespace Hts {

NlTaskManager::NlTaskManager(SliceMapping * mapping, RemoteState::MpiManager * commMgr,
                             NlPropagatorManager * propagatorMgr,
                             SeedInitializer * seedInitializer,
                             NlBasisUpdate * basisUpdateMgr,
                             PostProcessing::Manager * postProcessingMgr,
                             JumpConvergenceEvaluator * jumpCvgMgr, NonLinSeedDifferenceEvaluator::Manager * jumpEvaluatorMgr,
                             double projectorTolerance, IterationRank lastIteration) :
  TaskManager(IterationRank(0)),
  mapping_(mapping),
  commMgr_(commMgr),
  propagatorMgr_(propagatorMgr),
  seedInitializer_(seedInitializer),
  postProcessingMgr_(postProcessingMgr),
  jumpCvgMgr_(jumpCvgMgr),
  jumpEvaluatorMgr_(jumpEvaluatorMgr),
  projectorTolerance_(projectorTolerance),
  phase_(NULL),
  continuation_(&NlTaskManager::noop),
  localNetwork_(NULL),
  basisUpdateMgr_(basisUpdateMgr),
  projectionNetwork_(NULL),
  lastIteration_(lastIteration)
{
  initialize();
  scheduleInitialSeed();
}

void
NlTaskManager::iterationInc() {
  setIteration(iteration().next());
  if (iteration() < lastIteration_) {
    basisUpdateMgr_->globalUpdate()->mappingIs(*mapping_);
  }
  scheduleConvergence();
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
  NlDynamOps::Ptr dynOps = propagatorMgr_->dynamOpsNew();
  projectionNetwork_ = new NlProjectionNetwork(basisUpdateMgr_->globalUpdate(), dynOps.ptr(), commMgr_->vectorSize(), projectorTolerance_);

  // Parallel time-integration
  propagatorMgr_->postProcessingManagerIs(postProcessingMgr_.ptr());
  propagatorMgr_->concurrentBasisManagerIs(projectionNetwork_->concurrentMgr());
  propagatorMgr_->propagatedBasisManagerIs(projectionNetwork_->endBasisMgr());

  // Recurring tasks 
  localNetwork_ = new NlLocalNetwork(mapping_.ptr(),
                                     commMgr_.ptr(),
                                     propagatorMgr_.ptr(),
                                     projectionNetwork_->corrRedMgr(),
                                     projectionNetwork_->corrReconMgr(),
                                     projectionNetwork_->condensMgr(),
                                     projectionNetwork_->projBuildMgr(),
                                     jumpCvgMgr_.ptr(),
                                     jumpEvaluatorMgr_.ptr());
  localNetwork_->statusIs(LocalNetwork::ACTIVE);

  basisUpdateMgr_->globalUpdate()->seedGetterIs(localNetwork_->fullSeedGetter());
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
  setContinuation(&NlTaskManager::schedulePropagatedSeedSynchronization);
}

void
NlTaskManager::scheduleDataSharing() {
  TaskList dataSharing;
  if (iteration() < lastIteration_) {
    dataSharing.push_back(basisUpdateMgr_->globalUpdate());
    setPhase(phaseNew("Data Sharing", dataSharing));
    setContinuation(&NlTaskManager::enrichProjectionBasis);
  } else {
    projectionNetwork_->sharedProjectionBasis()->stateBasisDel();
    scheduleProjectionBasisCondensation();
  }
}

void
NlTaskManager::scheduleConvergence() {
  TaskList convergence;
  convergence.push_back(jumpCvgMgr_);
  setPhase(phaseNew("Convergence", convergence));
  setContinuation(&NlTaskManager::applyConvergence);
}

void
NlTaskManager::applyConvergence() {
  localNetwork_->applyConvergenceStatus();
  scheduleProjectionBuilding();
}

void
NlTaskManager::schedulePropagatedSeedSynchronization() {
  LocalNetwork::TaskList oldStyleList = localNetwork_->activePropagatedSeedSyncs();
  setPhase(phaseNew("Propagated Seed Synchronization", TaskList(oldStyleList.begin(), oldStyleList.end())));
  setContinuation(&NlTaskManager::scheduleJumpEvaluation);
}

void
NlTaskManager::scheduleJumpEvaluation() {
  LocalNetwork::TaskList oldStyleList = localNetwork_->activeJumpBuilders();
  setPhase(phaseNew("Jump Evaluation", TaskList(oldStyleList.begin(), oldStyleList.end())));
  setContinuation(&NlTaskManager::scheduleNothing);
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
  setContinuation(&NlTaskManager::scheduleDataSharing);
}

void
NlTaskManager::enrichProjectionBasis() {
  projectionNetwork_->sharedProjectionBasis()->firstStateBasisIs(basisUpdateMgr_->globalUpdate()->consolidatedBasis());
  scheduleProjectionBasisCondensation();
}

} /* end namespace Hts */ } /* end namespace Pita */
