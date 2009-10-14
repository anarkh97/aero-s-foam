#include "RemoteSeedInitializerServer.h"

#include <Comm.d/Communicator.h>

namespace Pita { namespace Hts {

RemoteSeedInitializerServer::RemoteSeedInitializerServer(Communicator * cc, SeedInitializer * si, SliceMapping * m) :
  clientCommunicator_(cc),
  baseInitializer_(si),
  mapping_(m),
  status_(READY),
  sBuffer_()
{}

void
RemoteSeedInitializerServer::statusIs(RemoteSeedInitializerServer::Status s) {
  if (status() == s)
    return;

  if (s == BUSY) {
    int seedCount = (mapping_->activeSlices().value() / 2) + 1;
    size_t stateSize = 2 * baseInitializer_->vectorSize();
    sBuffer_.sizeIs(stateSize * seedCount);
 
    for (int s = 0; s < seedCount; ++s) {
      DynamState seed = baseInitializer_->initialSeed(SliceRank(s));
      double * buffer = sBuffer_.array() + s * stateSize;
      bufferStateCopy(seed, buffer);
     
      CpuRank backwardTarget(-1);
      if (s != 0) {
        backwardTarget = mapping_->hostCpu(HalfSliceRank(s * 2 - 1));
        clientCommunicator_->sendTo(backwardTarget.value(), s, buffer, stateSize);
      }

      CpuRank forwardTarget(-1);
      if (s != seedCount - 1) {
        forwardTarget = mapping_->hostCpu(HalfSliceRank(s * 2));
        if (forwardTarget != backwardTarget) {
          clientCommunicator_->sendTo(forwardTarget.value(), s, buffer, stateSize);
        }
      }
    }

    clientCommunicator_->waitForAllReq();
  }
  
  status_ = READY;
}

} /* end namespace Hts */ } /* end namespace Hts */
