#include "RemoteSeedWriter.h"
#include "Comm.d/Communicator.h"
#include "Log.h"

namespace Pita {

template <>
void
RemoteSharedStateWriter<DynamState>::statusIs(Status s) {
  if (s != status()) {
    setStatus(s);
    if (s == BUSY) {
      if (originCpu().value() != timeComm_->myID()) {
        log() << "RemoteWriter for " << (target() ? target()->name() : "None") << " from Cpu " << originCpu()  << " is busy\n";
        size_t dim = 2 * vectorSize_ + 1;
        stateBuffer_.sizeIs(dim);
        int messageTag = originCpu().value(); // TODO messageTag
        timeComm_->recFrom(messageTag, stateBuffer_.array(), dim); // Blocking MPI_RECV
        log() << "RemoteWriter for " << (target() ? target()->name() : "None") << " from Cpu " << originCpu()  << " is done\n";
        if (target()) {
          target()->statusIs(Seed::Status(static_cast<int>(stateBuffer_[dim - 1])));
          target()->stateIs(DynamState(vectorSize_, stateBuffer_.array()));
        }
      }
      setStatus(READY);
    }
  }
}

} // namespace Pita
