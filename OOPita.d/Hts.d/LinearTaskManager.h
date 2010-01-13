#ifndef PITA_HTS_LINEARTASKMANAGER_H
#define PITA_HTS_LINEARTASKMANAGER_H

#include "../TaskManager.h"

#include "LocalNetwork.h"

#include "CorrectionNetworkImpl.h"
#include "../RemoteStateMpiImpl.h"

namespace Pita { namespace Hts {

class LinearTaskManager : public TaskManager {
public:
  EXPORT_PTRINTERFACE_TYPES(LinearTaskManager);

protected:
  LinearTaskManager(IterationRank initialIteration,
                    LocalNetwork * network,
                    CorrectionNetworkImpl * correctionMgr,
                    RemoteState::MpiManager * commMgr);

  LocalNetwork * network() { return network_.ptr(); }
  CorrectionNetworkImpl * correctionMgr() { return correctionMgr_.ptr(); }
  RemoteState::MpiManager * commMgr() { return commMgr_.ptr(); }

  void scheduleNormalIteration();
  void scheduleFinePropagation();
  void scheduleCorrection();

  void schedulePhase(const String & phaseName, const LocalNetwork::TaskList & networkTaskList);

private:
  LocalNetwork::Ptr network_;
  CorrectionNetworkImpl::Ptr correctionMgr_;
  RemoteState::MpiManager::Ptr commMgr_;
};

} /* end namespace Hts */ } /* end namespace Pita */

#endif /* PITA_HTS_LINEARTASKMANAGER_H */
