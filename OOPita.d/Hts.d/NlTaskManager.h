#ifndef PITA_HTS_NLTASKMANAGER_H
#define PITA_HTS_NLTASKMANAGER_H

#include "../TaskManager.h"

namespace Pita { namespace Hts {

class NlTaskManager : public TaskManager {
public:
  EXPORT_PTRINTERFACE_TYPES(NlTaskManager);

  virtual void iterationInc();
//protected:
  NlTaskManager();
protected:
  void scheduleIterationZero();
  void scheduleNormalIteration();
 
  void scheduleSeedInitialization(); 
  void scheduleFinePropagation();
  void scheduleCorrection();
};

} /* end namespace Hts */ } /* end namespace Pita */

#endif /* PITA_HTS_NLTASKMANAGER_H */
