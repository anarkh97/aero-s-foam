#include "NlTaskManager.h"

namespace Pita { namespace Hts {

NlTaskManager::NlTaskManager() :
  TaskManager(IterationRank(0))
{
  scheduleIterationZero();
}

void
NlTaskManager::iterationInc() {
  //phases().clear();
  
  IterationRank nextIteration = iteration().next();
  
  // TODO Pre-iteration 
  //convergedSlicesInc();
  scheduleNormalIteration();
  
  setIteration(nextIteration);
}

void
NlTaskManager::scheduleNormalIteration() {
  scheduleCorrection();
  scheduleFinePropagation();
}

void
NlTaskManager::scheduleIterationZero() {
  scheduleSeedInitialization();
  scheduleFinePropagation();
}

void
NlTaskManager::scheduleSeedInitialization() {
}

void
NlTaskManager::scheduleFinePropagation() {
}

void NlTaskManager::scheduleCorrection() {
}

} /* end namespace Hts */ } /* end namespace Pita */
