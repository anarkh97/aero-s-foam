#include "TaskManager.h"

#include <Utils.d/Connectivity.h>

namespace Pita {

TaskManager::TaskManager(TaskCount totTasks,  PoolCount availPools, TaskCount maxLoad) :
  totalTasks_(totTasks),
  maxWorkload_(maxLoad),
  availablePools_(availPools),
  completedTasks_(0),
  taskToPool_(NULL),
  poolToTasks_(NULL)
{
  int * arrayTasks = new int[totalTasks() + 1];
  for (int i = 0; i <= totalTasks(); ++i)
    arrayTasks[i] = i;

  int * targetPool = new int[totalTasks()];
  
  int ratio, remain;
  if (totalTasks() < availablePools() * maxWorkload()) {
    ratio = totalTasks() / availablePools();
    remain = totalTasks() - ratio * availablePools();
  } else {
    ratio = maxWorkload();
    remain = 0;
  }

  int currentPool = 0;
  int remainingSpotsOnCurrentPool = ratio + (remain > currentPool ? 1 : 0);
  for (int currentTask = 0; currentTask < totalTasks(); ++currentTask) {
    if (remainingSpotsOnCurrentPool <= 0) {
      currentPool = (currentPool + 1) % availablePools();
      remainingSpotsOnCurrentPool = ratio + (remain > currentPool ? 1 : 0);
    }
    targetPool[currentTask] = currentPool;
    --remainingSpotsOnCurrentPool;
  }

  taskToPool_ = new Connectivity(totalTasks(), arrayTasks, targetPool, 1);
  poolToTasks_ = taskToPool_->reverse();
  poolToTasks_->sortTargets(); 
  updateFirstWaitingTask();
}
  
TaskManager::~TaskManager() {
  delete taskToPool_;
  delete poolToTasks_;
}

void
TaskManager::completedTasksInc(TaskCount increment) {
  completedTasks_ += increment;
  updateFirstWaitingTask();
}

TaskCount
TaskManager::totalWorkload(PoolRank pr) const {
  return TaskCount(poolToTasks_->num(pr));
}

TaskCount
TaskManager::currentGlobalWorkload() const {
  return firstWaitingTask() - completedTasks();
}

TaskCount
TaskManager::currentWorkload(PoolRank pr) const {
  TaskCount result(0);
  for (TaskIterator it = tasks(pr); it && (*it < firstWaitingTask()); ++it) {
    if (*it >= completedTasks()) {
      ++result;
    }
  }
  return result;
}

TaskManager::TaskIterator
TaskManager::tasks(PoolRank pr) const {
  return TaskIterator(this, pr); 
}

PoolRank
TaskManager::pool(TaskRank tr) const {
  PoolRank result = (tr >= 0 && tr < totalTasks()) ?
                    taskToPool_->getTargetValue(tr) :
                    -1;
  return result;
}

TaskRank
TaskManager::firstCurrentTask() const { 
  return TaskRank(completedTasks());
}

void
TaskManager::updateFirstWaitingTask() {
  firstWaitingTask_ = std::min(totalTasks(), maxWorkload() *  availablePools() + completedTasks());
}

// TaskManager::TaskIterator

TaskManager::TaskIterator::TaskIterator(const TaskManager * parent, PoolRank pool) :
  parent_(parent),
  pool_(pool),
  offset_(0)
{
  if (parent_->totalWorkload(pool_) == TaskCount(0))
    parent_ = NULL;
}

TaskManager::TaskIterator &
TaskManager::TaskIterator::operator++() {
  if (parent_) {
    if ((++offset_) >= parent_->totalWorkload(pool_))
      parent_ = NULL;
  }
  return *this;
}

TaskManager::TaskIterator
TaskManager::TaskIterator::operator++(int) {
  TaskIterator temp(*this);
  ++(*this);
  return temp;
}

TaskRank
TaskManager::TaskIterator::operator*() const {
  TaskRank result = parent_ ?
                    *(parent_->poolToTasks_->operator[](pool_) + offset_) :
                    -1;
  return result;
}

TaskManager::TaskIterator::operator bool() const {
  return (parent_ != NULL);
}

// Output

OStream & operator<<(OStream & out, const TaskManager & tmgr) {
  for (int pool = 0; pool < tmgr.availablePools(); ++pool) {
    out << "# " << pool << " ->";
    for (TaskManager::TaskIterator it = tmgr.tasks(pool); it; ++it) {
      out << ' ' << *it;
    }
    out << '\n';
  }
 return out; 
}

} // end namespace Pita
