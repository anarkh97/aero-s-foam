#ifndef PITA_TASKMANAGER_H
#define PITA_TASKMANAGER_H

#include "Fwk.h"

class Connectivity;

namespace Pita {

typedef int TaskCount;
typedef int TaskRank;
typedef int PoolCount;
typedef int PoolRank;

class TaskManager : public Fwk::PtrInterface<TaskManager> {
public:
  EXPORT_PTRINTERFACE_TYPES(TaskManager);
  
  class TaskIterator;

  TaskCount totalTasks() const { return totalTasks_; }
  TaskCount maxWorkload() const { return maxWorkload_; }
  PoolCount availablePools() const { return availablePools_; }

  TaskCount completedTasks() const { return completedTasks_; }
  void completedTasksInc(TaskCount increment = TaskCount(1));

  TaskCount totalWorkload(PoolRank pr) const;
  TaskCount currentGlobalWorkload() const; 
  TaskCount currentWorkload(PoolRank pr) const; 

  TaskIterator tasks(PoolRank cr) const;
  PoolRank pool(TaskRank tr) const;

  TaskRank firstCurrentTask() const;
  TaskRank firstWaitingTask() const { return firstWaitingTask_; }

  static Ptr New(TaskCount totalTasks, PoolCount availablePools, TaskCount maxWorkload) {
    return new TaskManager(totalTasks, availablePools, maxWorkload);
  }

protected:
  TaskManager(TaskCount totalTasks, PoolCount availablePools, TaskCount maxWorkload);
  virtual ~TaskManager();

  void updateFirstWaitingTask();
  
  friend class TaskIterator;

private:
  TaskCount totalTasks_;
  TaskCount maxWorkload_;
  PoolCount availablePools_;
  TaskCount completedTasks_;

  Connectivity * taskToPool_;
  Connectivity * poolToTasks_;

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
  TaskIterator(const TaskManager * parent, PoolRank pool);

  friend class TaskManager;

private:
  TaskManager::PtrConst parent_;
  PoolRank pool_;
  TaskCount offset_;
};

OStream & operator<<(OStream & out, const TaskManager & tmgr);

} // end namespace Pita

#endif /* PITA_TASKMANAGER_H */
