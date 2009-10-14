#ifndef PITA_HTS_REMOTESEEDINITIALIZERSERVER_H
#define PITA_HTS_REMOTESEEDINITIALIZERSERVER_H

#include "../SeedInitializer.h"
#include "SliceMapping.h"
#include "../SimpleBuffer.h"

class Communicator;

namespace Pita { namespace Hts {

class RemoteSeedInitializerServer : public Fwk::PtrInterface<RemoteSeedInitializerServer> {
public:
  EXPORT_PTRINTERFACE_TYPES(RemoteSeedInitializerServer);

  enum Status {
    READY = 0,
    BUSY
  };

  Status status() const { return status_; }
  void statusIs(Status s);

  static Ptr New(Communicator * clientCommunicator, SeedInitializer * baseInitializer, SliceMapping * mapping) {
    return new RemoteSeedInitializerServer(clientCommunicator, baseInitializer, mapping);
  }

protected:
  RemoteSeedInitializerServer(Communicator * cc, SeedInitializer * si, SliceMapping * m);

private:
  Communicator * clientCommunicator_;
  SeedInitializer::Ptr baseInitializer_;
  SliceMapping::Ptr mapping_;

  Status status_;

  SimpleBuffer<double> sBuffer_;
};

} // end namespace Hts
} // end namespace Pita

#endif /* PITA_HTS_REMOTESEEDINITIALIZERSERVER_H */
