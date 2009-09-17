#ifndef PITA_TASKMANAGER_H
#define PITA_TASKMANAGER_H

#include "Fwk.h"

class Connectivity;

namespace Pita {

typedef int TaskCount;
typedef int TaskRank;
typedef int WorkerCount;
typedef int WorkerRank;

class TaskManager : public Fwk::PtrInterface<TaskManager> {
public:
  EXPORT_PTRINTERFACE_TYPES(TaskManager);
  
  class TaskIterator;

  TaskCount totalTasks() const { return totalTasks_; }
  TaskCount maxWorkload() const { return maxWorkload_; }
  WorkerCount availableWorkers() const { return availableWorkers_; }

  TaskCount completedTasks() const { return completedTasks_; }
  void completedTasksInc(TaskCount increment = TaskCount(1));

  TaskCount totalWorkload(WorkerRank pr) const;
  TaskCount currentGlobalWorkload() const; 
  TaskCount currentWorkload(WorkerRank pr) const; 

  TaskIterator tasks(WorkerRank cr) const;
  WorkerRank worker(TaskRank tr) const;

  TaskRank firstCurrentTask() const;
  TaskRank firstWaitingTask() const { return firstWaitingTask_; }

  static Ptr New(TaskCount totalTasks, WorkerCount availableWorkers, TaskCount maxWorkload) {
    return new TaskManager(totalTasks, availableWorkers, maxWorkload);
  }

protected:
  TaskManager(TaskCount totalTasks, WorkerCount availableWorkers, TaskCount maxWorkload);
  virtual ~TaskManager();

  void updateFirstWaitingTask();
  
  friend class TaskIterator;

private:
  TaskCount totalTasks_;
  TaskCount maxWorkload_;
  WorkerCount availableWorkers_;
  TaskCount completedTasks_;

  Connectivity * taskToWorker_;
  Connectivity * workerToTasks_;

  TaskRank firstWaitingTask_;

  DISALLOW_COPY_AND_ASSIGN(TaskManager);
};

class TaskManager::TaskIterator {
public:
  TaskIterator & operator++();
  TaskIterator operator++(int);
  TaskRank operator*() const;
  operator bool() const;

  // Use default copy constructor, assignement operator

protected:
  TaskIterator(const TaskManager * parent, WorkerRank worker);

  friend class TaskManager;

private:
  TaskManager::PtrConst parent_;
  WorkerRank worker_;
  TaskCount offset_;
};

OStream & operator<<(OStream & out, const TaskManager & tmgr);

} // end namespace Pita

#endif /* PITA_TASKMANAGER_H */
