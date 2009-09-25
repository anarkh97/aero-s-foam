#ifndef PITA_HTS_REMOTECOARSECORRECTIONSERVER_H
#define PITA_HTS_REMOTECOARSECORRECTIONSERVER_H

#include "Fwk.h"
#include "Types.h"
#include "SliceMapping.h"
#include "../RemoteDynamPropagatorServer.h"

#include "../Activity.h"
#include <map>

namespace Pita { namespace Hts {

class RemoteCoarseCorrectionServer : public Fwk::PtrInterface<RemoteCoarseCorrectionServer> {
public:
  EXPORT_PTRINTERFACE_TYPES(RemoteCoarseCorrectionServer);

  enum Status {
    IDLE = 0,
    ACTIVE
  };

  Status status() const { return status_; }
  RemoteDynamPropagatorServer * server() const { return server_.ptr(); }
  const SliceMapping * mapping() const { return mapping_.ptr(); }
  PhaseRank correctionPhase() const { return correctionPhase_; }

  void statusIs(Status s);

  RemoteCoarseCorrectionServer(
      RemoteDynamPropagatorServer * server,
      const SliceMapping * mapping,
      PhaseRank correctionPhase);

protected:
  class CorrectionReactor : public Activity::Notifiee {
  public:
    EXPORT_PTRINTERFACE_TYPES(CorrectionReactor);

    virtual void onStatus();

    CorrectionReactor(Activity * notifier, RemoteCoarseCorrectionServer * parent);
  private:
    RemoteCoarseCorrectionServer * parent_;
  };

private:
  Status status_;
  RemoteDynamPropagatorServer::Ptr server_;
  SliceMapping::PtrConst mapping_;
  PhaseRank correctionPhase_;

  CorrectionReactor::Ptr correctionReactor_;

  DISALLOW_COPY_AND_ASSIGN(RemoteCoarseCorrectionServer);
};

} /* end namespace Hts */ } /* end namespace Pita */

#endif /* PITA_HTS_REMOTECOARSECORRECTIONSERVER_H */
