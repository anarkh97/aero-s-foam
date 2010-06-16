#include "NonHomogeneousTaskManager.h"

namespace Pita { namespace Hts {

class SimpleInitialSeed : public NamedTask {
public:
  EXPORT_PTRINTERFACE_TYPES(SimpleInitialSeed);

  virtual void iterationIs(IterationRank i) {
    target_->stateIs(state_);
    target_->statusIs(status_);
    target_->iterationIs(i);
  }

  SimpleInitialSeed(Seed * target, DynamState state, Seed::Status status) :
    NamedTask("Seed initialization " + toString(target->name())),
    target_(target),
    state_(state),
    status_(status)
  {}

private:
  Seed::Ptr target_;
  DynamState state_;
  Seed::Status status_;
};

NonHomogeneousTaskManager::NonHomogeneousTaskManager(LinearLocalNetwork * network,
                                                     LinearProjectionNetworkImpl * correctionMgr,
                                                     RemoteState::MpiManager * commMgr,
                                                     DynamState initialCondition) :
  LinearTaskManager(IterationRank(-1), network, correctionMgr, commMgr),
  initialCondition_(initialCondition)
{
  schedulePreIteration();
  updatePhaseIt();
}

void
NonHomogeneousTaskManager::iterationInc() {
  phases().clear();
  
  IterationRank nextIteration = iteration().next();

  if (nextIteration == IterationRank(0)) {
    network()->convergedSlicesInc();
    scheduleIterationZero();
  } else {
    correctionMgr()->prepareProjection();
    network()->convergedSlicesInc();
    scheduleNormalIteration();
  }

  updatePhaseIt();  
  setIteration(nextIteration);
}

void
NonHomogeneousTaskManager::schedulePreIteration() {
  scheduleBasicSeedInitialization();
  schedulePhase("Affine term precomputation", network()->halfTimeSlices());
}

void
NonHomogeneousTaskManager::scheduleBasicSeedInitialization() {
  TaskList taskList;

  DynamState zeroState(initialCondition_.vectorSize(), 0.0);

  LinearLocalNetwork::SeedMap mainSeeds = network()->mainSeeds();
  for (LinearLocalNetwork::SeedMap::iterator it = mainSeeds.begin();
      it != mainSeeds.end();
      ++it) {
  
    Seed * target = it->second.ptr(); 
    DynamState state = zeroState;
    Seed::Status status = (it->first.value() % 2 == 0) ? Seed::ACTIVE : Seed::SPECIAL;

    if (target->name() == "M0")  {
      state = initialCondition_;
      status = Seed::CONVERGED;
    }

    taskList.push_back(new SimpleInitialSeed(target, state, status));
  }

  Phase::Ptr zeroSeedInitialization = phaseNew("Basic seed initialization", taskList);
  phases().push_back(zeroSeedInitialization);
}

void
NonHomogeneousTaskManager::scheduleIterationZero() {
  schedulePhase("Propagated seed synchronization", network()->activeLeftSeedSyncs());
  schedulePhase("Jumps", network()->activeJumpProjectors());
  schedulePhase("Coarse propagation", network()->activeCoarseTimeSlices());
  schedulePhase("Correction synchronization", network()->activeFullCorrectionSyncs());
  schedulePhase("Seed update", network()->activeSeedAssemblers());
  scheduleFinePropagation(); 
}

} /* end namespace Hts */ } /* end namespace Pta */
