#ifndef PITA_HTS_NONHOMOGENEOUSTASKMANAGER_H
#define PITA_HTS_NONHOMOGENEOUSTASKMANAGER_H

#include "LinearTaskManager.h"

namespace Pita { namespace Hts {

class NonHomogeneousTaskManager : public LinearTaskManager {
public:
  EXPORT_PTRINTERFACE_TYPES(NonHomogeneousTaskManager);

  virtual void iterationInc(); //overriden

  NonHomogeneousTaskManager(LinearLocalNetwork * network,
                            LinearProjectionNetworkImpl * correctionMgr,
                            RemoteState::MpiManager * commMgr,
                            DynamState initialCondition);

protected:
  void schedulePreIteration();
  void scheduleIterationZero();
  void scheduleBasicSeedInitialization();

private:
  DynamState initialCondition_;
};

} /* end namespace Hts */ } /* end namespace Pta */

#endif /* PITA_HTS_NONHOMOGENEOUSTASKMANAGER_H */
