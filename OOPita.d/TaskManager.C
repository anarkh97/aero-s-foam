#include "TaskManager.h"

#include <Utils.d/Connectivity.h>

namespace Pita {

TaskManager::TaskManager(TaskCount totTasks, WorkerCount availWorkers, TaskCount maxLoad) :
  totalTasks_(totTasks),
  maxWorkload_(maxLoad),
  availableWorkers_(availWorkers),
  completedTasks_(0),
  taskToWorker_(NULL),
  workerToTasks_(NULL)
{
  int * arrayTasks = new int[totalTasks() + 1];
  for (int i = 0; i <= totalTasks(); ++i)
    arrayTasks[i] = i;

  int * targetWorker = new int[totalTasks()];
  
  int ratio, remain;
  if (totalTasks() < availableWorkers() * maxWorkload()) {
    ratio = totalTasks() / availableWorkers();
    remain = totalTasks() - ratio * availableWorkers();
  } else {
    ratio = maxWorkload();
    remain = 0;
  }

  int currentWorker = 0;
  int remainingSpotsForCurrentWorker = ratio + (remain > currentWorker ? 1 : 0);
  for (int currentTask = 0; currentTask < totalTasks(); ++currentTask) {
    if (remainingSpotsForCurrentWorker <= 0) {
      currentWorker = (currentWorker + 1) % availableWorkers();
      remainingSpotsForCurrentWorker = ratio + (remain > currentWorker ? 1 : 0);
    }
    targetWorker[currentTask] = currentWorker;
    --remainingSpotsForCurrentWorker;
  }

  taskToWorker_ = new Connectivity(totalTasks(), arrayTasks, targetWorker, 1);
  workerToTasks_ = taskToWorker_->reverse();
  workerToTasks_->sortTargets(); 
  updateFirstWaitingTask();
}
  
TaskManager::~TaskManager() {
  delete taskToWorker_;
  delete workerToTasks_;
}

void
TaskManager::completedTasksInc(TaskCount increment) {
  completedTasks_ += increment;
  updateFirstWaitingTask();
}

TaskCount
TaskManager::totalWorkload(WorkerRank pr) const {
  return TaskCount(workerToTasks_->num(pr));
}

TaskCount
TaskManager::currentGlobalWorkload() const {
  return firstWaitingTask() - completedTasks();
}

TaskCount
TaskManager::currentWorkload(WorkerRank pr) const {
  TaskCount result(0);
  for (TaskIterator it = tasks(pr); it && (*it < firstWaitingTask()); ++it) {
    if (*it >= completedTasks()) {
      ++result;
    }
  }
  return result;
}

TaskManager::TaskIterator
TaskManager::tasks(WorkerRank pr) const {
  return TaskIterator(this, pr); 
}

WorkerRank
TaskManager::worker(TaskRank tr) const {
  return (tr >= 0 && tr < totalTasks())    ?
         taskToWorker_->getTargetValue(tr) :
         -1;
}

TaskRank
TaskManager::firstCurrentTask() const { 
  return TaskRank(completedTasks());
}

void
TaskManager::updateFirstWaitingTask() {
  firstWaitingTask_ = std::min(totalTasks(), maxWorkload() *  availableWorkers() + completedTasks());
}

// TaskManager::TaskIterator

TaskManager::TaskIterator::TaskIterator(const TaskManager * parent, WorkerRank worker) :
  parent_(parent),
  worker_(worker),
  offset_(0)
{
  if (parent_->totalWorkload(worker_) == TaskCount(0))
    parent_ = NULL;
}

TaskManager::TaskIterator &
TaskManager::TaskIterator::operator++() {
  if (parent_) {
    if ((++offset_) >= parent_->totalWorkload(worker_))
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
                    *(parent_->workerToTasks_->operator[](worker_) + offset_) :
                    -1;
  return result;
}

TaskManager::TaskIterator::operator bool() const {
  return (parent_ != NULL);
}

// Output

OStream & operator<<(OStream & out, const TaskManager & tmgr) {
  for (int worker = 0; worker < tmgr.availableWorkers(); ++worker) {
    out << "# " << worker << " ->";
    for (TaskManager::TaskIterator it = tmgr.tasks(worker); it; ++it) {
      out << ' ' << *it;
    }
    out << '\n';
  }
 return out; 
}

} // end namespace Pita
