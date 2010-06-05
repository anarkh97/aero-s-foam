#ifndef PITA_HTS_REMOTECOARSESERVER_H
#define PITA_HTS_REMOTECOARSESERVER_H

#include "SliceMapping.h"

class Communicator;

namespace Pita { namespace Hts {

class RemoteCoarseServer : public Fwk::PtrInterface<RemoteCoarseServer> {
public:
  EXPORT_PTRINTERFACE_TYPES(RemoteCoarseServer);

  enum Status {
    READY = 0,
    BUSY
  };

  Status status() const { return status_; }
  virtual void statusIs(Status s) = 0;

  Communicator * clientCommunicator() const { return clientCommunicator_; }
  const SliceMapping * clientMapping() const { return clientMapping_.ptr(); }

protected:
  RemoteCoarseServer(Communicator * clientComm, const SliceMapping * clientMapping) :
    clientCommunicator_(clientComm),
    clientMapping_(clientMapping),
    status_(READY)
  {}

  void setStatus(Status s) { status_ = s; }

private:
  Communicator * clientCommunicator_;
  SliceMapping::PtrConst clientMapping_;

  Status status_;
};

} /* end namespace Hts */ } /* end namespace Pita */

#endif /* PITA_HTS_REMOTECOARSESERVER_H */
