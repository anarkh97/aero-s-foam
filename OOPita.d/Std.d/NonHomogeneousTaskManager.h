#ifndef PITA_STD_NONHOMOGENEOUSTASKMANAGER_H
#define PITA_STD_NONHOMOGENEOUSTASKMANAGER_H

#include "../TaskManager.h"

#include "SliceMapping.h"

#include "LinearPropagatorManager.h"
#include "../SeedInitializer.h"

#include "LinearProjectionNetwork.h"

#include "../JumpProjector.h"
#include "../CorrectionPropagator.h"
#include "../UpdatedSeedAssembler.h"

#include "../SeedDifferenceEvaluator.h"

#include "../Seed.h"

#include "../RemoteStateMpiImpl.h"

#include <map>

namespace Pita { namespace Std {

class NonHomogeneousTaskManager : public TaskManager {
public:
  EXPORT_PTRINTERFACE_TYPES(NonHomogeneousTaskManager);

  virtual Phase * phase(); // overriden
  virtual void phaseInc(); // overriden
  virtual void iterationInc(); // overriden

  NonHomogeneousTaskManager(SliceMapping * mapping, CpuRank localCpu, RemoteState::MpiManager * commMgr,
                            DynamState initialState, LinearPropagatorManager * propagatorMgr, LinearProjectionNetwork * projectionMgr,
                            JumpProjector::Manager * jumpProjMgr, CorrectionPropagator<Vector>::Manager * corrPropMgr,
                            CorrectionPropagator<DynamState>::Manager * fullCorrPropMgr, UpdatedSeedAssembler::Manager * seedUpMgr);

protected:
  // Initialization
  void initialize();
  void addLocalSlice(SliceRank slice);
  void addPrecedingRemoteSlice(SliceRank slice);
  void addFollowingRemoteSlice(SliceRank slice);

  // Execution
  void setPhase(Phase * p) { phase_ = p; }

  typedef void (NonHomogeneousTaskManager::*Continuation)();
  void setContinuation(Continuation c) { continuation_ = c; }

  void noop() {}
  void scheduleNothing();
  
  void scheduleInitialization();
  void schedulePrecomputation();
  void scheduleAffineTermSynchronization();
  void scheduleJumpBuilding();
  void scheduleCoarseCorrectionPropagation();

  void scheduleFinePropagation();
  void scheduleProjectionBuilding();
  void schedulePropagatedSeedSynchronization();
  void scheduleCorrectionPropagation();
  void scheduleJumpProjection();
  void scheduleSeedUpdate();

  void checkConvergence();

private:
  SliceMapping::Ptr mapping_;
  CpuRank localCpu_;

  SeedInitializer::Ptr initializer_;
  LinearPropagatorManager::Ptr propagatorMgr_;
  
  LinearProjectionNetwork::Ptr projectionMgr_;
  
  JumpProjector::Manager::Ptr jumpProjMgr_;
  UpdatedSeedAssembler::Manager::Ptr seedUpMgr_;
  CorrectionPropagator<Vector>::Manager::Ptr corrPropMgr_;
  
  JumpBuilder::Manager::Ptr jumpBuildMgr_;
  CorrectionPropagator<DynamState>::Manager::Ptr fullCorrPropMgr_;

  Seed::Manager::Ptr seedMgr_;
  ReducedSeed::Manager::Ptr redSeedMgr_;

  RemoteState::MpiManager::Ptr commMgr_;

  struct SliceTasks {
    NamedTask::Ptr propagatedSeedSynchronization;
    NamedTask::Ptr jumpProjection;
    NamedTask::Ptr correctionPropagation;
    NamedTask::Ptr seedUpdate;
    NamedTask::Ptr finePropagation;
  };
  typedef NamedTask::Ptr SliceTasks::*SliceTaskItem;

  typedef std::map<SliceRank, SliceTasks> SliceTaskMap;
  SliceTaskMap recurrentTasks_;

  void fillTaskList(SliceTaskItem item, TaskList & target);

  LinSeedDifferenceEvaluator::Manager::Ptr jumpEvalMgr_;

  Phase::Ptr phase_;

  Continuation continuation_;
  
  DISALLOW_COPY_AND_ASSIGN(NonHomogeneousTaskManager);
};

} /* end namespace Std */ } /* end namespace Pita */

#endif /* PITA_STD_NONHOMOGENEOUSTASKMANAGER_H */
