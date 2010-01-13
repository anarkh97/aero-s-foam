#ifndef PITA_TASKMANAGER_H
#define PITA_TASKMANAGER_H

#include "Fwk.h"
#include "Types.h"

#include "NamedTask.h"

#include <vector>

namespace Pita {

class TaskManager : public Fwk::PtrInterface<TaskManager> {
public:
  EXPORT_PTRINTERFACE_TYPES(TaskManager);

  // Iteration control
  IterationRank iteration() const { return iteration_; }
  virtual void iterationInc() = 0;

  // Phases
  class Phase;

protected:
  // Auxiliary type
  template <typename E>
  class Iterator {
  public:
    E * operator*() { return it_->ptr(); }
    E * operator->() { return it_->ptr(); }
    Iterator<E> & operator++() { ++it_; return *this; }
    Iterator<E> operator++(int) { Iterator<E> temp(*this); this->operator++(); return temp; }
    operator bool() const { return it_ != it_end_; } 

    explicit Iterator(std::vector<Fwk::Ptr<E> > & container) :
      it_(container.begin()),
      it_end_(container.end())
    {}

  private:
    typename std::vector<Fwk::Ptr<E> >::iterator it_;
    typename std::vector<Fwk::Ptr<E> >::iterator it_end_;
  };


public:
  class Phase : public NamedInterface {
  public:
    EXPORT_PTRINTERFACE_TYPES(Phase);

    typedef Iterator<NamedTask> TaskIterator;

    TaskIterator task() { return TaskIterator(task_); }

  protected:
    Phase(const String & name, const std::vector<NamedTask::Ptr> & taskList) :
      NamedInterface(name),
      task_(taskList)
    {}

    friend class TaskManager;

  private:
    typedef std::vector<NamedTask::Ptr> TaskList;
    TaskList task_; 

    DISALLOW_COPY_AND_ASSIGN(Phase);
  };
  
  typedef Iterator<Phase> PhaseIterator;

  PhaseIterator phase() { return PhaseIterator(phase_); }

protected:
  explicit TaskManager(IterationRank initialIter) :
    iteration_(initialIter)
  {}

  void setIteration(IterationRank i) { iteration_ = i; }

  typedef std::vector<Fwk::Ptr<Phase> > PhaseList;
  typedef Phase::TaskList TaskList;

  static Phase * phaseNew(const String & name, const TaskList & taskList) {
    return new Phase(name, taskList);
  }

  PhaseList & phases() { return phase_; }

private:
  IterationRank iteration_;
  PhaseList phase_;

  DISALLOW_COPY_AND_ASSIGN(TaskManager);
};


} /* end namespace Pita */

#endif /* PITA_TASKMANAGER_H */
