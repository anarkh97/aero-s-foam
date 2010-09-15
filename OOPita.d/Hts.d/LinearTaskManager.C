#include "LinearTaskManager.h"

namespace Pita { namespace Hts {

class ProjectionBasis : public NamedTask {
public:
  void iterationIs(IterationRank i) {
    correctionMgr_->buildProjection();
    size_t reducedBasisSize = correctionMgr_->reducedBasisSize();
    commMgr_->reducedStateSizeIs(reducedBasisSize);
  }

  ProjectionBasis(LinearProjectionNetworkImpl * correctionMgr, RemoteState::MpiManager * commMgr) :
    NamedTask("Projection building"),
    correctionMgr_(correctionMgr),
    commMgr_(commMgr)
  {}

private:
  LinearProjectionNetworkImpl::Ptr correctionMgr_;
  RemoteState::MpiManager::Ptr commMgr_;
};

LinearTaskManager::LinearTaskManager(IterationRank initialIteration,
                                     LinearLocalNetwork * network,
                                     JumpConvergenceEvaluator * jumpCvgMgr,
                                     LinearProjectionNetworkImpl * correctionMgr,
                                     RemoteState::MpiManager * commMgr) :
  TaskManager(initialIteration),
  network_(network),
  jumpCvgMgr_(jumpCvgMgr),
  correctionMgr_(correctionMgr),
  commMgr_(commMgr),
  phase_(),

  phaseIt_(new HtsPhaseIteratorImpl(phase_))
{}

/*class ApplyConvergence : public NamedTask {
public:
  EXPORT_PTRINTERFACE_TYPES(ApplyConvergence);
  virtual void iterationIs(IterationRank iter); // overriden

  explicit ApplyConvergence(LinearLocalNetwork * network) :
    NamedTask("Apply convergence"),
    network_(network)
  {}

private:
  LinearLocalNetwork::Ptr network_;
};

void
ApplyConvergence::iterationIs(IterationRank iter) {
  network_->convergedSlicesInc(); // TODO bad naming
  setIteration(iter);
}*/

void
LinearTaskManager::scheduleNormalIteration() {
  log() << "Executing: " << jumpCvgMgr()->name() << "\n";
  jumpCvgMgr()->iterationIs(iteration().next()); // TODO BETTER
  network()->applyConvergenceStatus();

  scheduleCorrection();
  scheduleFinePropagation();
}

void
LinearTaskManager::scheduleFinePropagation() {
  schedulePhase("Fine propagation", network()->activeHalfTimeSlices());
  schedulePhase("Propagated seed synchronization", network()->activeLeftSeedSyncs());
  schedulePhase("Jump evaluation", network()->activeJumpAssemblers());
}

//void
//LinearTaskManager::schedulePreCorrection() {
  // Convergence
  /*{
    TaskList tasks;
    tasks.push_back(jumpCvgMgr());
    tasks.push_back(new ApplyConvergence(network()));
    Phase::Ptr convergence = phaseNew("Convergence", tasks);
    phases().push_back(convergence);
  }*/
//}

void
LinearTaskManager::scheduleCorrection() {
  // Projection Basis
  {
    TaskList projectionBuilding;
    projectionBuilding.push_back(new ProjectionBasis(correctionMgr_.ptr(), commMgr_.ptr()));
    Phase::Ptr projectionBasis = phaseNew("Projection basis", projectionBuilding);
    phases().push_back(projectionBasis);
  }
  schedulePhase("Jump Projection", network()->activeJumpProjectors());
  schedulePhase("Correction", network()->activeFullTimeSlices());
  schedulePhase("Correction synchronization", network()->activeCorrectionSyncs());
  schedulePhase("Seed update", network()->activeSeedAssemblers());
}

void
LinearTaskManager::schedulePhase(const String & phaseName, const LinearLocalNetwork::TaskList & networkTaskList) {
  TaskList taskList(networkTaskList.begin(), networkTaskList.end());

  Phase::Ptr projectionBasis = phaseNew(phaseName, taskList);
  phases().push_back(projectionBasis);
}

} /* end namespace Hts */ } /* end namespace Pita */
