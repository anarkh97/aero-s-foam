#ifndef PITA_TIMEDEXECUTION_H
#define PITA_TIMEDEXECUTION_H

#include "Fwk.h"

#include "TaskManager.h"

namespace Pita {

class TimedExecution : public Fwk::PtrInterface<TimedExecution> {
public:
  EXPORT_PTRINTERFACE_TYPES(TimedExecution);

  IterationRank currentIteration() const { return currentIteration_; }
  void targetIterationIs(IterationRank ti);

  TaskManager * taskManager() const { return taskMgr_.ptr(); }

  static Ptr New(TaskManager * taskMgr) {
    return new TimedExecution(taskMgr);
  }

protected:
  explicit TimedExecution(TaskManager * taskMgr);

private:
  TaskManager::Ptr taskMgr_;
  IterationRank currentIteration_;

  DISALLOW_COPY_AND_ASSIGN(TimedExecution);
};

} /* end namespace Pita */

#endif /* PITA_TIMEDEXECUTION_H */
