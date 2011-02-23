#ifndef PITA_STD_LINEARTASKMANAGER_H
#define PITA_STD_LINEARTASKMANAGER_H

#include "../TaskManager.h"

#include "SliceMapping.h"

#include "LinearPropagatorManager.h"
#include "LinearProjectionNetwork.h"
#include "JumpConvergenceEvaluator.h"

#include "../JumpBuilder.h"
#include "../JumpProjection.h"
#include "../ReducedCorrectionPropagatorImpl.h"
#include "../UpdatedSeedAssemblerImpl.h"

#include "../RemoteStateMpiImpl.h"
#include "../SeedDifferenceEvaluator.h"

#include "../SeedInitializer.h"

#include "../InitialSeedTask.h"
#include "../IncrementalPropagation.h"
#include "../RemoteStateTask.h"

#include <map>

namespace Pita { namespace Std {

// Note: Impl must be a subtype of LinearTaskManager<Impl>
template <typename Impl>
class LinearTaskManager : public TaskManager {
public:
  EXPORT_PTRINTERFACE_TYPES(LinearTaskManager);

  virtual Phase * phase(); // overriden
  virtual void phaseInc(); // overriden
  virtual void iterationInc(); // overriden

  const SliceMapping * mapping() const { return mapping_.ptr(); }
  SliceMapping * mapping() { return mapping_.ptr(); }

  CpuRank localCpu() const { return localCpu_; }

protected:
  LinearTaskManager(IterationRank initialIteration,
                    SliceMapping * mapping,
                    RemoteState::MpiManager * commMgr,
                    LinearPropagatorManager * propagatorMgr,
                    LinearProjectionNetwork * projectionMgr,
                    JumpConvergenceEvaluator * jumpCvgEval,
                    LinSeedDifferenceEvaluator::Manager * jumpOutMgr,
                    SeedInitializer * initializer);

  // Initialization
  void initialize();
  void addLocalSlice(SliceRank slice);
  void addPrecedingRemoteSlice(SliceRank slice);
  void addFollowingRemoteSlice(SliceRank slice);

  // Execution
  void setPhase(Phase * p) { phase_ = p; }

  typedef void (Impl::*Continuation)();
  void setContinuation(Continuation c) { continuation_ = c; }

  void noop() {}
  void scheduleNothing();
  
  void scheduleInitialization();
  void scheduleFinePropagation();
  void schedulePropagatedSeedSynchronization();
  void scheduleJumpBuilding();
  void scheduleProjectionBuilding();
  void scheduleConvergence();
  void scheduleJumpProjection();
  void scheduleCorrectionPropagation();
  void scheduleSeedUpdate();

  void applyConvergence();

private:
  SliceMapping::Ptr mapping_;
  CpuRank localCpu_;

  LinearPropagatorManager::Ptr propagatorMgr_;
  LinearProjectionNetwork::Ptr projectionMgr_;
  JumpConvergenceEvaluator::Ptr jumpCvgEval_;
 
  JumpBuilder::Manager::Ptr jumpBuildMgr_; 
  JumpProjection::Manager::Ptr jumpProjMgr_;
  UpdatedSeedAssembler::Manager::Ptr seedUpMgr_;
  CorrectionPropagator<Vector>::Manager::Ptr corrPropMgr_;

protected:
  Seed::Manager::Ptr seedMgr_;
  ReducedSeed::Manager::Ptr redSeedMgr_;
  RemoteState::MpiManager::Ptr commMgr_;

private:
  LinSeedDifferenceEvaluator::Manager::Ptr jumpOutMgr_;
  
  SeedInitializer::Ptr initializer_;

  Phase::Ptr phase_;
  Continuation continuation_;
 
protected: 
  struct SliceTasks {
    NamedTask::Ptr propagatedSeedSynchronization;
    NamedTask::Ptr jumpBuilding;
    NamedTask::Ptr jumpProjection;
    NamedTask::Ptr correctionPropagation;
    NamedTask::Ptr seedUpdate;
    NamedTask::Ptr finePropagation;
  };
  typedef NamedTask::Ptr SliceTasks::*SliceTaskItem;

  void setRecurrentPhase(const String & name, SliceTaskItem item);

  typedef std::map<SliceRank, SliceTasks> SliceTaskMap;
  SliceTaskMap recurrentTasks_;

  void fillTaskList(SliceTaskItem item, TaskList & target);

  DISALLOW_COPY_AND_ASSIGN(LinearTaskManager);
};

template <typename Impl>
LinearTaskManager<Impl>::LinearTaskManager(IterationRank initialIteration,
                                           SliceMapping * mapping,
                                           RemoteState::MpiManager * commMgr,
                                           LinearPropagatorManager * propagatorMgr,
                                           LinearProjectionNetwork * projectionMgr,
                                           JumpConvergenceEvaluator * jumpCvgEval,
                                           LinSeedDifferenceEvaluator::Manager * jumpOutMgr,
                                           SeedInitializer * initializer) :
  TaskManager(initialIteration),
  mapping_(mapping),
  localCpu_(commMgr->localCpu()),
  initializer_(initializer),
  propagatorMgr_(propagatorMgr),
  projectionMgr_(projectionMgr),
  jumpBuildMgr_(JumpBuilder::ManagerImpl::New()),
  jumpProjMgr_(JumpProjection::Manager::New(projectionMgr->projectionBasis())),
  seedUpMgr_(UpdatedSeedAssemblerImpl::Manager::New(projectionMgr->propagatedBasis())),
  corrPropMgr_(ReducedCorrectionPropagatorImpl::Manager::New(projectionMgr->reprojectionMatrix(), projectionMgr->normalMatrixSolver())),
  seedMgr_(Seed::Manager::New()),
  redSeedMgr_(ReducedSeed::Manager::New()),
  commMgr_(commMgr),
  jumpCvgEval_(jumpCvgEval),
  jumpOutMgr_(jumpOutMgr), 
  phase_(NULL),
  continuation_(&Impl::noop)
{
  initialize();
  scheduleInitialization();
}

template <typename Impl>
void
LinearTaskManager<Impl>::iterationInc() {
  scheduleProjectionBuilding();
  setIteration(iteration().next());
}

template <typename Impl>
TaskManager::Phase *
LinearTaskManager<Impl>::phase() {
  return phase_.ptr();
}

template <typename Impl>
void
LinearTaskManager<Impl>::phaseInc() {
  Impl * self = static_cast<Impl*>(this); // Safe if Impl is a subtype of LinearTaskManager<Impl>
  (self->*continuation_)();
}

template <typename Impl>
void
LinearTaskManager<Impl>::scheduleNothing() {
  setPhase(NULL);
  setContinuation(&Impl::noop);
}

template <typename Impl>
void
LinearTaskManager<Impl>::scheduleInitialization() {
  TaskList initialSeedInformation;

  for (SliceMapping::SliceIterator sl_it = mapping_->hostedSlice(localCpu_); sl_it; ++sl_it) {
    SliceRank slice = *sl_it;
    Seed::Ptr targetSeed = seedMgr_->instance(toString(SeedId(MAIN_SEED, slice)));
    Seed::Status initialSeedStatus = (slice == SliceRank(0)) ? Seed::CONVERGED : Seed::ACTIVE;
    NamedTask::Ptr task = new InitialSeedTask(targetSeed.ptr(), initializer_.ptr(), slice, initialSeedStatus);
    initialSeedInformation.push_back(task);
  }
  
  setPhase(phaseNew("Initial Seed Information", initialSeedInformation));
}

template <typename Impl>
inline
void
LinearTaskManager<Impl>::setRecurrentPhase(const String & name, SliceTaskItem item) {
  TaskList tasksInPhase;
  fillTaskList(item, tasksInPhase);
  setPhase(phaseNew(name, tasksInPhase));
}

template <typename Impl>
void
LinearTaskManager<Impl>::scheduleFinePropagation() {
  setRecurrentPhase("Fine Propagation", &SliceTasks::finePropagation);
  setContinuation(&Impl::schedulePropagatedSeedSynchronization);
}

template <typename Impl>
void
LinearTaskManager<Impl>::schedulePropagatedSeedSynchronization() {
  setRecurrentPhase("Propagated Seed Synchronization", &SliceTasks::propagatedSeedSynchronization);
  setContinuation(&Impl::scheduleJumpBuilding);
}

template <typename Impl>
void
LinearTaskManager<Impl>::scheduleJumpBuilding() {
  setRecurrentPhase("Jump Evaluation", &SliceTasks::jumpBuilding);
  setContinuation(&Impl::scheduleNothing);
}

template <typename Impl>
void
LinearTaskManager<Impl>::scheduleProjectionBuilding() {
  TaskList projBuilding(1, projectionMgr_->projectionTaskNew());
  setPhase(phaseNew("Build Projection", projBuilding));
  setContinuation(&Impl::scheduleConvergence);
}

template <typename Impl>
void
LinearTaskManager<Impl>::scheduleConvergence() {
  TaskList checkCvg(1, jumpCvgEval_); 
  setPhase(phaseNew("Check Convergence", checkCvg));
  setContinuation(&Impl::scheduleJumpProjection);
}

template <typename Impl>
void
LinearTaskManager<Impl>::scheduleJumpProjection() {
  applyConvergence();

  setRecurrentPhase("Jump Projection", &SliceTasks::jumpProjection);
  setContinuation(&Impl::scheduleCorrectionPropagation);
}

template <typename Impl>
void
LinearTaskManager<Impl>::scheduleCorrectionPropagation() {
  commMgr_->reducedStateSizeIs(projectionMgr_->reducedBasisSize());

  setRecurrentPhase("Correction Propagation", &SliceTasks::correctionPropagation);
  setContinuation(&Impl::scheduleSeedUpdate);
}

template <typename Impl>
void
LinearTaskManager<Impl>::scheduleSeedUpdate() {
  setRecurrentPhase("Seed Update", &SliceTasks::seedUpdate);
  setContinuation(&Impl::scheduleFinePropagation);
}

template <typename Impl>
void
LinearTaskManager<Impl>::applyConvergence() {
  SliceRank firstActiveSlice = mapping_->firstActiveSlice();
  recurrentTasks_.erase(recurrentTasks_.begin(), recurrentTasks_.lower_bound(firstActiveSlice));
 
  String seedStr = toString(SeedId(SEED_CORRECTION, firstActiveSlice)); 
  
  Seed::Ptr fullLeadingCorrection = seedMgr_->instance(seedStr);
  if (fullLeadingCorrection) {
    fullLeadingCorrection->statusIs(Seed::INACTIVE);
  }
  
  ReducedSeed::Ptr reducedLeadingCorrection = redSeedMgr_->instance(seedStr);
  if (reducedLeadingCorrection) {
    reducedLeadingCorrection->statusIs(Seed::INACTIVE);
  }
  
  if (localCpu_ == mapping_->hostCpu(firstActiveSlice.previous()) && localCpu_ != mapping_->hostCpu(firstActiveSlice)) {
    recurrentTasks_[firstActiveSlice].correctionPropagation = NULL;
  }
}

template <typename Impl>
void
LinearTaskManager<Impl>::initialize() {
  commMgr_->vectorSizeIs(initializer_->vectorSize());

  for (SliceMapping::SliceIterator sl_it = mapping_->hostedSlice(localCpu_); sl_it; ++sl_it) {
    SliceRank slice = *sl_it;
    SliceRank previousSlice = slice.previous();
    SliceRank nextSlice = slice.next();

    if (previousSlice >= mapping_->firstActiveSlice() && mapping_->hostCpu(previousSlice) != localCpu_) {
      addPrecedingRemoteSlice(previousSlice);
    }

    addLocalSlice(slice);

    if (nextSlice < mapping_->firstInactiveSlice() && mapping_->hostCpu(nextSlice) != localCpu_) {
      addFollowingRemoteSlice(nextSlice);
    }
  }
}

template <typename Impl>
void
LinearTaskManager<Impl>::addLocalSlice(SliceRank slice) {
  String sliceStr = toString(slice);

  // Fine propagation
  {
    Seed::Ptr seed = seedMgr_->instanceNew(toString(SeedId(MAIN_SEED, slice))); 
    Seed::Ptr propSeed = seedMgr_->instanceNew(toString(SeedId(PROPAGATED_SEED, slice.next())));

    AffineDynamPropagator::Ptr propagator = propagatorMgr_->instanceNew(slice);
    IncrementalPropagation::Ptr task = IncrementalPropagation::New(String("Propagate ") + sliceStr, propagator.ptr());
    task->seedIs(seed.ptr());
    task->propagatedSeedIs(propSeed.ptr());

    recurrentTasks_[slice].finePropagation = task;
  }

  // Jump building, convergence and output
  if (slice > SliceRank(0)) {
    Seed::Ptr seed = seedMgr_->instance(toString(SeedId(MAIN_SEED, slice))); 
    Seed::Ptr prevPropSeed = seedMgr_->instance(toString(SeedId(PROPAGATED_SEED, slice)));
    Seed::Ptr jump = seedMgr_->instanceNew(toString(SeedId(SEED_JUMP, slice)));

    JumpBuilder::Ptr task = jumpBuildMgr_->instanceNew(sliceStr);

    task->predictedSeedIs(seed.ptr());
    task->actualSeedIs(prevPropSeed.ptr());
    task->seedJumpIs(jump.ptr());

    recurrentTasks_[slice.previous()].jumpBuilding = task;
  
    jumpCvgEval_->localJumpIs(slice, jump.ptr());

    if (jumpOutMgr_) { 
      Seed::Ptr propagatedSeed = seedMgr_->instance(toString(SeedId(PROPAGATED_SEED, slice)));
      LinSeedDifferenceEvaluator::Ptr eval = jumpOutMgr_->instanceNew(jump.ptr());
      eval->referenceSeedIs(propagatedSeed.ptr());
    }
  }
 
  // Jump projection
  {
    Seed::Ptr jump = seedMgr_->instance(toString(SeedId(SEED_JUMP, slice)));
    ReducedSeed::Ptr redJump = redSeedMgr_->instanceNew(toString(SeedId(SEED_JUMP, slice)));

    JumpProjection::Ptr task = jumpProjMgr_->instanceNew(sliceStr);

    task->seedJumpIs(jump.ptr());
    task->reducedSeedJumpIs(redJump.ptr());

    recurrentTasks_[slice].jumpProjection = task;
  }

  // Correction propagation
  {
    ReducedSeed::Ptr redJump = redSeedMgr_->instance(toString(SeedId(SEED_JUMP, slice)));
    ReducedSeed::Ptr correction = redSeedMgr_->instance(toString(SeedId(SEED_CORRECTION, slice)));
    ReducedSeed::Ptr nextCorrection = redSeedMgr_->instanceNew(toString(SeedId(SEED_CORRECTION, slice.next())));

    CorrectionPropagator<Vector>::Ptr task = corrPropMgr_->instanceNew(sliceStr);

    task->jumpIs(redJump.ptr());
    task->correctionIs(correction.ptr());
    task->nextCorrectionIs(nextCorrection.ptr());

    recurrentTasks_[slice].correctionPropagation = task;
  }

  // Seed update
  {
    Seed::Ptr seed = seedMgr_->instance(toString(SeedId(MAIN_SEED, slice)));
    Seed::Ptr fullCorrection = seedMgr_->instanceNew(toString(SeedId(SEED_CORRECTION, slice)));
    Seed::Ptr prevPropSeed = seedMgr_->instance(toString(SeedId(PROPAGATED_SEED, slice)));
    ReducedSeed::Ptr reducedCorrection = redSeedMgr_->instance(toString(SeedId(SEED_CORRECTION, slice)));

    UpdatedSeedAssembler::Ptr task = seedUpMgr_->instanceNew(sliceStr);

    task->updatedSeedIs(seed.ptr());
    task->propagatedSeedIs(prevPropSeed.ptr());
    task->correctionIs(fullCorrection.ptr());
    task->correctionComponentsIs(reducedCorrection.ptr());

    recurrentTasks_[slice].seedUpdate = task;
  }
}

template <typename Impl>
void
LinearTaskManager<Impl>::addPrecedingRemoteSlice(SliceRank slice) {
  CpuRank peer = mapping_->hostCpu(slice);
  
  // 1) Receive Prop seed sync from previous
  {
    Seed::Ptr propagatedSeed = seedMgr_->instanceNew(toString(SeedId(PROPAGATED_SEED, slice.next())));
    RemoteState::Writer<DynamState>::Ptr activity = commMgr_->writerNew(propagatedSeed.ptr(), peer);
    RemoteStateTask::Ptr task = RemoteStateTask::New("Receive " + propagatedSeed->name(), activity.ptr());
    recurrentTasks_[slice].propagatedSeedSynchronization = task;
  }

  // 2) Receive Correction from previous
  {
    ReducedSeed::Ptr reducedCorrection = redSeedMgr_->instanceNew(toString(SeedId(SEED_CORRECTION, slice.next())));
    RemoteState::Writer<Vector>::Ptr activity = commMgr_->writerNew(reducedCorrection.ptr(), peer);
    RemoteStateTask::Ptr task = RemoteStateTask::New("Receive " + reducedCorrection->name(), activity.ptr());
    recurrentTasks_[slice].correctionPropagation = task;
  }
}

template <typename Impl>
void
LinearTaskManager<Impl>::addFollowingRemoteSlice(SliceRank slice) {
  CpuRank peer = mapping_->hostCpu(slice);
  
  // 1) Send Prop seed sync to next
  {
    Seed::Ptr propagatedSeed = seedMgr_->instance(toString(SeedId(PROPAGATED_SEED, slice)));
    RemoteState::Reader<DynamState>::Ptr activity = commMgr_->readerNew(propagatedSeed.ptr(), peer);
    RemoteStateTask::Ptr task = RemoteStateTask::New("Send " + propagatedSeed->name(), activity.ptr());
    recurrentTasks_[slice.previous()].propagatedSeedSynchronization = task;
  }

  // 2) Send Correction to next
  {
    ReducedSeed::Ptr reducedCorrection = redSeedMgr_->instance(toString(SeedId(SEED_CORRECTION, slice)));
    RemoteState::Reader<Vector>::Ptr activity = commMgr_->readerNew(reducedCorrection.ptr(), peer);
    RemoteStateTask::Ptr task = RemoteStateTask::New("Send " + reducedCorrection->name(), activity.ptr());
    recurrentTasks_[slice].correctionPropagation = task;
  }
}

template <typename Impl>
void
LinearTaskManager<Impl>::fillTaskList(SliceTaskItem item, TaskList & target) {
  typename SliceTaskMap::const_iterator it_end = recurrentTasks_.end();
  for (typename SliceTaskMap::const_iterator it = recurrentTasks_.begin(); it != it_end; ++it) {
    NamedTask::Ptr task = (it->second).*item;
    if (task) {
      target.push_back(task);
    }
  }
}

} /* end namespace Std */ } /* end namespace Pita */

#endif /* PITA_STD_LINEARTASKMANAGER_H */
