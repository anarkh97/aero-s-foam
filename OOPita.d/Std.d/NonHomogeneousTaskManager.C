#include "NonHomogeneousTaskManager.h"

#include "../SimpleSeedInitializer.h"

#include "../InitialSeedTask.h"
#include "IncrementalPropagation.h"
#include "../RemoteStateTask.h"

namespace Pita { namespace Std {

NonHomogeneousTaskManager::NonHomogeneousTaskManager(SliceMapping * mapping,
                                                     CpuRank localCpu,
                                                     RemoteState::MpiManager * commMgr,
                                                     DynamState initialState,
                                                     LinearPropagatorManager * propagatorMgr,
                                                     LinearProjectionNetwork * projectionMgr,
                                                     JumpProjection::Manager * jumpProjMgr,
                                                     CorrectionPropagator<Vector>::Manager * corrPropMgr,
                                                     CorrectionPropagator<DynamState>::Manager * fullCorrPropMgr,
                                                     UpdatedSeedAssembler::Manager * seedUpMgr,
                                                     JumpConvergenceEvaluator * jumpCvgEval,
                                                     LinSeedDifferenceEvaluator::Manager * jumpOutMgr) :
  TaskManager(IterationRank(-1)),
  mapping_(mapping),
  localCpu_(localCpu),
  initializer_(SimpleSeedInitializer::New(initialState)),
  propagatorMgr_(propagatorMgr),
  projectionMgr_(projectionMgr),
  jumpProjMgr_(jumpProjMgr),
  seedUpMgr_(seedUpMgr),
  corrPropMgr_(corrPropMgr),
  jumpBuildMgr_(JumpBuilder::ManagerImpl::New()),
  fullCorrPropMgr_(fullCorrPropMgr),
  seedMgr_(Seed::Manager::New()),
  redSeedMgr_(ReducedSeed::Manager::New()),
  commMgr_(commMgr),
  jumpCvgEval_(jumpCvgEval),
  jumpOutMgr_(jumpOutMgr),
  phase_(NULL),
  continuation_(&NonHomogeneousTaskManager::noop)
{
  initialize();
  scheduleInitialization();
}

void
NonHomogeneousTaskManager::iterationInc() {
  setIteration(iteration().next());
  if (iteration() == IterationRank(0)) {
    scheduleCoarseCorrectionPropagation();
  } else {
    scheduleProjectionBuilding();
  }
}


NonHomogeneousTaskManager::Phase *
NonHomogeneousTaskManager::phase() {
  return phase_.ptr();
}

void
NonHomogeneousTaskManager::phaseInc() {
  (this->*continuation_)();
}

void
NonHomogeneousTaskManager::scheduleNothing() {
  setPhase(NULL);
  setContinuation(&NonHomogeneousTaskManager::noop);
}

void
NonHomogeneousTaskManager::scheduleInitialization() {
  TaskList initialSeedInformation;

  for (SliceMapping::SliceIterator sl_it = mapping_->hostedSlice(localCpu_); sl_it; ++sl_it) {
    SliceRank slice = *sl_it;
    Seed::Ptr targetSeed = seedMgr_->instance(toString(SeedId(MAIN_SEED, slice)));
    Seed::Status initialSeedStatus = (slice == SliceRank(0)) ? Seed::CONVERGED : Seed::ACTIVE;
    NamedTask::Ptr task = new InitialSeedTask(targetSeed.ptr(), initializer_.ptr(), slice, initialSeedStatus);
    initialSeedInformation.push_back(task);
  }
  
  setPhase(phaseNew("Initial Seed Information", initialSeedInformation));
  setContinuation(&NonHomogeneousTaskManager::schedulePrecomputation);
}

void
NonHomogeneousTaskManager::schedulePrecomputation() {
  TaskList finePropagation;
  fillTaskList(&SliceTasks::finePropagation, finePropagation);
  setPhase(phaseNew("Affine Term Precomputation", finePropagation));

  setContinuation(&NonHomogeneousTaskManager::scheduleAffineTermSynchronization);
}

void
NonHomogeneousTaskManager::scheduleAffineTermSynchronization() {
  TaskList propagatedSeedSynchronization;
  fillTaskList(&SliceTasks::propagatedSeedSynchronization, propagatedSeedSynchronization);
  setPhase(phaseNew("Affine Term Synchronization", propagatedSeedSynchronization));

  setContinuation(&NonHomogeneousTaskManager::scheduleTrivialJumpBuilding);
}

void
NonHomogeneousTaskManager::scheduleTrivialJumpBuilding() {
  TaskList jumpBuilding;
  fillTaskList(&SliceTasks::jumpBuilding, jumpBuilding);
  setPhase(phaseNew("Trivial Jump Evaluation", jumpBuilding));

  setContinuation(&NonHomogeneousTaskManager::scheduleNothing);
}

void
NonHomogeneousTaskManager::scheduleCoarseCorrectionPropagation() {
  mapping_->convergedSlicesInc();
  applyConvergence();
  
  TaskList correctionPropagation;
 
  for (SliceMapping::SliceIterator sl_it = mapping_->hostedSlice(localCpu_); sl_it; ++sl_it) {
    SliceRank slice = *sl_it;
    SliceRank previousSlice = slice.previous();
    SliceRank nextSlice = slice.next();

    if (slice < mapping_->firstActiveSlice()) {
      continue;
    }

    if (slice >= mapping_->firstInactiveSlice()) {
      break;
    }

    // Receive previous  
    if (previousSlice >= mapping_->firstActiveSlice() && mapping_->hostCpu(previousSlice) != localCpu_) {
      CpuRank peer = mapping_->hostCpu(previousSlice);
      
      Seed::Ptr correction = seedMgr_->instance(toString(SeedId(SEED_CORRECTION, slice)));
      RemoteState::Writer<DynamState>::Ptr activity = commMgr_->writerNew(correction.ptr(), peer);
      RemoteStateTask::Ptr task = RemoteStateTask::New("Receive " + correction->name(), activity.ptr());
      
      correctionPropagation.push_back(task);
    } else {
      seedMgr_->instance(toString(SeedId(SEED_CORRECTION, slice)))->statusIs(Seed::INACTIVE);
    }

    // Local propagation
    {
      Seed::Ptr jump = seedMgr_->instance(toString(SeedId(SEED_JUMP, slice)));
      Seed::Ptr correction = seedMgr_->instance(toString(SeedId(SEED_CORRECTION, slice)));
      Seed::Ptr nextCorrection = seedMgr_->instance(toString(SeedId(SEED_CORRECTION, nextSlice)));
      if (!nextCorrection) {
        nextCorrection = seedMgr_->instanceNew(toString(SeedId(SEED_CORRECTION, nextSlice)));
      }

      CorrectionPropagator<DynamState>::Ptr task = fullCorrPropMgr_->instanceNew(toString(slice));

      task->jumpIs(jump.ptr());
      task->correctionIs(correction.ptr());
      task->nextCorrectionIs(nextCorrection.ptr());

      correctionPropagation.push_back(task);
    }

    // Send next
    if (nextSlice < mapping_->firstInactiveSlice() && mapping_->hostCpu(nextSlice) != localCpu_) {
      CpuRank peer = mapping_->hostCpu(nextSlice);

      Seed::Ptr correction = seedMgr_->instance(toString(SeedId(SEED_CORRECTION, nextSlice)));
      RemoteState::Reader<DynamState>::Ptr activity = commMgr_->readerNew(correction.ptr(), peer);
      RemoteStateTask::Ptr task = RemoteStateTask::New("Send " + correction->name(), activity.ptr());
      
      correctionPropagation.push_back(task);
    }
  }
  
  setPhase(phaseNew("Coarse Correction Propagation", correctionPropagation));
  
  setContinuation(&NonHomogeneousTaskManager::scheduleSeedUpdate);
}

void
NonHomogeneousTaskManager::scheduleFinePropagation() {
  TaskList finePropagation;
  fillTaskList(&SliceTasks::finePropagation, finePropagation);
  setPhase(phaseNew("Fine Propagation", finePropagation));

  setContinuation(&NonHomogeneousTaskManager::schedulePropagatedSeedSynchronization);
}

void
NonHomogeneousTaskManager::schedulePropagatedSeedSynchronization() {
  TaskList propagatedSeedSynchronization;
  fillTaskList(&SliceTasks::propagatedSeedSynchronization, propagatedSeedSynchronization);
  setPhase(phaseNew("Propagated Seed Synchronization", propagatedSeedSynchronization));

  setContinuation(&NonHomogeneousTaskManager::scheduleJumpBuilding);
}

void
NonHomogeneousTaskManager::scheduleJumpBuilding() {
  TaskList jumpBuilding;
  fillTaskList(&SliceTasks::jumpBuilding, jumpBuilding);
  setPhase(phaseNew("Jump Evaluation", jumpBuilding));

  setContinuation(&NonHomogeneousTaskManager::scheduleNothing);
}

void
NonHomogeneousTaskManager::scheduleProjectionBuilding() {
  TaskList projBuilding; 
  projBuilding.push_back(projectionMgr_->projectionTaskNew());
  setPhase(phaseNew("Build Projection", projBuilding));
 
  setContinuation(&NonHomogeneousTaskManager::scheduleConvergence);
}

void
NonHomogeneousTaskManager::scheduleConvergence() {
  TaskList checkCvg; 
  checkCvg.push_back(jumpCvgEval_.ptr());
  setPhase(phaseNew("Check Convergence", checkCvg));
 
  setContinuation(&NonHomogeneousTaskManager::scheduleJumpProjection);
}

void
NonHomogeneousTaskManager::scheduleJumpProjection() {
  applyConvergence();

  TaskList jumpProjection;
  fillTaskList(&SliceTasks::jumpProjection, jumpProjection);
  setPhase(phaseNew("Jump Projection", jumpProjection));

  setContinuation(&NonHomogeneousTaskManager::scheduleCorrectionPropagation);
}

void
NonHomogeneousTaskManager::scheduleCorrectionPropagation() {
  commMgr_->reducedStateSizeIs(projectionMgr_->reducedBasisSize());
  TaskList correctionPropagation;
  fillTaskList(&SliceTasks::correctionPropagation, correctionPropagation);
  setPhase(phaseNew("Correction Propagation", correctionPropagation));
  
  setContinuation(&NonHomogeneousTaskManager::scheduleSeedUpdate);
}

void
NonHomogeneousTaskManager::scheduleSeedUpdate() {
  TaskList seedUpdate;
  fillTaskList(&SliceTasks::seedUpdate, seedUpdate);
  setPhase(phaseNew("Seed Update", seedUpdate));
  
  setContinuation(&NonHomogeneousTaskManager::scheduleFinePropagation);
}

void
NonHomogeneousTaskManager::applyConvergence() {
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

void
NonHomogeneousTaskManager::initialize() {
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

void
NonHomogeneousTaskManager::addLocalSlice(SliceRank slice) {
  String sliceStr = toString(slice);

  // Fine propagation
  {
    Seed::Ptr seed = seedMgr_->instanceNew(toString(SeedId(MAIN_SEED, slice))); 
    Seed::Ptr propSeed = seedMgr_->instanceNew(toString(SeedId(PROPAGATED_SEED, slice.next())));

    AffineDynamPropagator::Ptr propagator = propagatorMgr_->instanceNew(slice);
    propagator->constantTermIs(AffineDynamPropagator::NONHOMOGENEOUS);
    IncrementalPropagation::Ptr task = IncrementalPropagation::New(String("Fine Propagation ") + sliceStr, propagator.ptr());
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

void
NonHomogeneousTaskManager::addPrecedingRemoteSlice(SliceRank slice) {
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

void
NonHomogeneousTaskManager::addFollowingRemoteSlice(SliceRank slice) {
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

void
NonHomogeneousTaskManager::fillTaskList(SliceTaskItem item, TaskList & target) {
  SliceTaskMap::const_iterator it_end = recurrentTasks_.end();
  for (SliceTaskMap::const_iterator it = recurrentTasks_.begin(); it != it_end; ++it) {
    NamedTask::Ptr task = (it->second).*item;
    if (task) {
      target.push_back(task);
    }
  }
}


} /* end namespace Std */ } /* end namespace Pita */
