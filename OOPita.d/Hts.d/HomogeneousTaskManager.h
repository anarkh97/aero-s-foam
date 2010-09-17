#ifndef PITA_HTS_HOMOGENEOUSTASKMANAGER_H
#define PITA_HTS_HOMOGENEOUSTASKMANAGER_H

#include "LinearTaskManager.h"

#include "../SeedInitializer.h"

namespace Pita { namespace Hts {

class HomogeneousTaskManager : public LinearTaskManager {
public:
  EXPORT_PTRINTERFACE_TYPES(HomogeneousTaskManager);

  // overriden
  virtual void iterationInc();

  explicit HomogeneousTaskManager(LinearLocalNetwork * network,
                                  SeedInitializer * initializer,
                                  JumpConvergenceEvaluator * jumpCvgMgr,
                                  LinearProjectionNetworkImpl * correctionMgr,
                                  RemoteState::MpiManager * commMgr);

protected:
  void scheduleIterationZero();
  void scheduleSeedInitialization();

private:
  SeedInitializer::Ptr initializer_;
};

} /* end namespace Hts */ } /* end namespace Pita */

#endif /* PITA_HTS_HOMOGENEOUSTASKMANAGER_H */
