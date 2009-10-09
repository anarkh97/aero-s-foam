#ifndef PITA_HTS_REMOTECOARSECORRECTIONSERVER_H
#define PITA_HTS_REMOTECOARSECORRECTIONSERVER_H

#include "Fwk.h"
#include "Types.h"
#include "SliceMapping.h"
#include "../RemoteDynamPropagatorServer.h"

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

  void statusIs(Status s);

  RemoteCoarseCorrectionServer(
      RemoteDynamPropagatorServer * server,
      const SliceMapping * mapping,
      PhaseRank correctionPhase);

private:
  Status status_;
  RemoteDynamPropagatorServer::Ptr server_;
  SliceMapping::PtrConst mapping_;

  DISALLOW_COPY_AND_ASSIGN(RemoteCoarseCorrectionServer);
};

} /* end namespace Hts */ } /* end namespace Pita */

#endif /* PITA_HTS_REMOTECOARSECORRECTIONSERVER_H */
