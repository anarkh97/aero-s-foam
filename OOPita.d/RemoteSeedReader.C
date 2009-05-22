#include "RemoteSeedReader.h"

// MPI Communication
#include <Comm.d/Communicator.h>

#include "Log.h"

namespace Pita {

RemoteSeedReader::RemoteSeedReader(Communicator * comm) :
  Seed::NotifieeConst(NULL),
  targetCpu_(),
  status_(READY),
  timeComm_(comm),
  stateBuffer_()
{}

void
RemoteSeedReader::onState() {
  this->statusIs(BUSY);
}

void
RemoteSeedReader::statusIs(Status s) {
  if (s != status()) {
    setStatus(s);
    if (status() == BUSY) {
      if (timeComm_->myID() != targetCpu().value() && targetCpu() != CpuRank(-1)) {
        DynamState state = notifier()->state();
        size_t dim = 2 * state.vectorSize() + 1;
        stateBuffer_.sizeIs(dim);
        bufferStateCopy(state, stateBuffer_.array());
        stateBuffer_[dim - 1] = static_cast<double>(notifier()->status()); 
        int messageTag = timeComm_->myID(); // TODO determine meaningful messageId
        log() << "RemoteReader for " << (notifier() ? notifier()->name() : "None") << " to Cpu " << targetCpu()  << " is busy\n";
        timeComm_->sendTo(targetCpu().value(), messageTag, stateBuffer_.array(), dim);
        timeComm_->waitForAllReq(); // TODO delay waitForAllReq
        log() << "RemoteReader for " << (notifier() ? notifier()->name() : "None") << " to Cpu " << targetCpu()  << " is done\n";
        // TODO
      }
      setStatus(READY);
    }
  }
}

} // end namespace Pita
