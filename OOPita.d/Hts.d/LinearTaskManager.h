#ifndef PITA_HTS_LINEARTASKMANAGER_H
#define PITA_HTS_LINEARTASKMANAGER_H

#include "../TaskManager.h"

#include "LinearLocalNetwork.h"

#include "LinearProjectionNetworkImpl.h"
#include "../RemoteStateMpiImpl.h"

namespace Pita { namespace Hts {

class LinearTaskManager : public TaskManager {
public:
  EXPORT_PTRINTERFACE_TYPES(LinearTaskManager);

protected:
  LinearTaskManager(IterationRank initialIteration,
                    LinearLocalNetwork * network,
                    LinearProjectionNetworkImpl * correctionMgr,
                    RemoteState::MpiManager * commMgr);

  LinearLocalNetwork * network() { return network_.ptr(); }
  LinearProjectionNetworkImpl * correctionMgr() { return correctionMgr_.ptr(); }
  RemoteState::MpiManager * commMgr() { return commMgr_.ptr(); }

  void scheduleNormalIteration();
  void scheduleFinePropagation();
  void scheduleCorrection();

  void schedulePhase(const String & phaseName, const LinearLocalNetwork::TaskList & networkTaskList);

private:
  LinearLocalNetwork::Ptr network_;
  LinearProjectionNetworkImpl::Ptr correctionMgr_;
  RemoteState::MpiManager::Ptr commMgr_;
};

} /* end namespace Hts */ } /* end namespace Pita */

#endif /* PITA_HTS_LINEARTASKMANAGER_H */
