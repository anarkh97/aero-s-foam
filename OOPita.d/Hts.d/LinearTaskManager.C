#include "LinearTaskManager.h"

namespace Pita { namespace Hts {

class ProjectionBasis : public NamedTask {
public:
  void iterationIs(IterationRank i) {
    correctionMgr_->buildProjection();
    size_t reducedBasisSize = correctionMgr_->reducedBasisSize();
    log() << "*** ReducedBasisSize = " << reducedBasisSize << "\n";
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
                                     LinearProjectionNetworkImpl * correctionMgr,
                                     RemoteState::MpiManager * commMgr) :
  TaskManager(initialIteration),
  network_(network),
  correctionMgr_(correctionMgr),
  commMgr_(commMgr)
{}

void
LinearTaskManager::scheduleNormalIteration() {
  scheduleCorrection();
  scheduleFinePropagation();
}

void
LinearTaskManager::scheduleFinePropagation() {
  schedulePhase("Fine propagation", network()->activeHalfTimeSlices());
}

void
LinearTaskManager::scheduleCorrection() {
  // Projection Basis
  TaskList taskList;
  taskList.push_back(new ProjectionBasis(correctionMgr_.ptr(), commMgr_.ptr()));
  Phase::Ptr projectionBasis = phaseNew("Projection basis", taskList);
  phases().push_back(projectionBasis);

  schedulePhase("Propagated seed synchronization", network()->activeLeftSeedSyncs());
  schedulePhase("Jumps", network()->activeJumpProjectors());
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
