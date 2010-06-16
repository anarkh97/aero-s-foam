#include "TimedExecution.h"

namespace Pita {

TimedExecution::TimedExecution(TaskManager * taskMgr) :
  taskMgr_(taskMgr),
  currentIteration_(taskMgr->iteration())
{}

void
TimedExecution::targetIterationIs(IterationRank lastIteration) {
  while (currentIteration() <= lastIteration) {
    log() << "-> Iteration: " << currentIteration() << "\n";

    while (TaskManager::Phase::Ptr currentPhase = taskMgr_->phase()) {
      log() << "  -> Phase: " << currentPhase->name() << "\n";
      for (TaskManager::Phase::TaskIterator task_it = currentPhase->task(); task_it; ++task_it) {
        NamedTask::Ptr currentTask = *task_it;
        log() << "    -> Task: " << currentTask->name() << "\n";
        currentTask->iterationIs(currentIteration());
      } 
      taskMgr_->phaseInc();
    }

    taskMgr_->iterationInc();
    currentIteration_ = taskMgr_->iteration();
  }
}

} /* end namespace Pita */
