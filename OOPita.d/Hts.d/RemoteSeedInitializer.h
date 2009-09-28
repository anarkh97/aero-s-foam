#ifndef PITA_HTS_REMOTESEEDINITIALIZER_H
#define PITA_HTS_REMOTESEEDINITIALIZER_H

#include "../SeedInitializer.h"
#include "SliceMapping.h"
#include "../SimpleBuffer.h"

class Communicator;

namespace Pita { namespace Hts {

class RemoteSeedInitializer : public Fwk::PtrInterface<RemoteSeedInitializer> {
public:
  EXPORT_PTRINTERFACE_TYPES(RemoteSeedInitializer);

  enum Status {
    READY = 0,
    BUSY
  };

  Status status() const { return status_; }
  void statusIs(Status s);

  static Ptr New(Communicator * clientCommunicator, SeedInitializer * baseInitializer, SliceMapping * mapping) {
    return new RemoteSeedInitializer(clientCommunicator, baseInitializer, mapping);
  }

protected:
  RemoteSeedInitializer(Communicator * cc, SeedInitializer * si, SliceMapping * m);

private:
  Communicator * clientCommunicator_;
  SeedInitializer::Ptr baseInitializer_;
  SliceMapping::Ptr mapping_;

  Status status_;

  SimpleBuffer<double> sBuffer_;
};

} // end namespace Hts
} // end namespace Pita

#endif /* PITA_HTS_REMOTESEEDINITIALIZER_H */
