#ifndef PITA_HTS_NONHOMOGENEOUSTASKMANAGER_H
#define PITA_HTS_NONHOMOGENEOUSTASKMANAGER_H

#include "LinearTaskManager.h"
#include "../SeedInitializer.h"

namespace Pita { namespace Hts {

class NonHomogeneousTaskManager : public LinearTaskManager {
public:
  EXPORT_PTRINTERFACE_TYPES(NonHomogeneousTaskManager);

  virtual void iterationInc(); //overriden

  NonHomogeneousTaskManager(LinearLocalNetwork * network,
                            SeedInitializer * seedInit,
                            JumpConvergenceEvaluator * jumpCvgMgr,
                            LinearProjectionNetworkImpl * correctionMgr,
                            RemoteState::MpiManager * commMgr);

protected:
  void schedulePreIteration();
  void scheduleIterationZero();
  void scheduleBasicSeedInitialization();

private:
  SeedInitializer::Ptr seedInit_;
};

} /* end namespace Hts */ } /* end namespace Pta */

#endif /* PITA_HTS_NONHOMOGENEOUSTASKMANAGER_H */
