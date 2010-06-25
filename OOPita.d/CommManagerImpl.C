#include "CommManagerImpl.h"
#include "DynamStateBasisWrapper.h"

// MPI Communication
#include <Comm.d/Communicator.h>

namespace Pita {

CommManagerImpl::CommManagerImpl(TimeSliceMapping * sliceMap, CpuRank myCpuRank, CpuRank nextCpuRank, CpuRank previousCpuRank, Communicator * timeComm) :
  CommManager(sliceMap, myCpuRank, nextCpuRank, previousCpuRank),
  timeCommunicator_(timeComm)
{}

CommManagerImpl::Ptr
CommManagerImpl::New(TimeSliceMapping * sliceMap, Communicator * timeComm) {
  CpuCount availableCpus(sliceMap->availableCpus());
  CpuRank myCpuRank(timeComm->myID());
  return new CommManagerImpl(sliceMap,
                             myCpuRank,
                             CpuRank((myCpuRank.value() + 1) % availableCpus.value()),
                             CpuRank((myCpuRank.value() + availableCpus.value() - 1) % availableCpus.value()),
                             timeComm);
}

void
CommManagerImpl::localBroadcastBasisIs(DynamStateBasis::Ptr basis) {
  setLocalBroadcastBasis(basis.ptr());
  
  // Setup data gBuffer
  size_t stateSize = 2 * basis->vectorSize();
  size_t localStateCount = basis->stateCount();
  size_t totalStateCount = sliceMap()->activeSlices().value();
  gBuffer_.sizeIs(stateSize * totalStateCount);

  // Setup parameters for global MPI communication
  int numCpus = availableCpus().value(); 
  mpiParameters_.sizeIs(2 * numCpus);
  int * recv_counts = mpiParameters_.array();
  int * displacements = mpiParameters_.array() + numCpus;

  //log() << "recv_counts: ";
  for (int i = 0; i < numCpus; ++i) {
    recv_counts[i] = stateSize * sliceMap()->activeSlicesOnCpu(CpuRank(i)).value();
    //log() << recv_counts[i] << " ";
  }
  //log() << "\n";
 
  //log() << "displacements: 0 "; 
  displacements[0] = 0;
  for (int i = 1; i < numCpus; ++i) {
    displacements[i] = displacements[i-1] + recv_counts[i-1];
    //log() << displacements[i] << " ";
  }
  //log() << "\n";

  // Fill the buffer
  unsigned int myCpu = localCpuRank().value();
  size_t sliceShift = displacements[myCpu];
  for (size_t i = 0; i < localStateCount; ++i) {
    bufferStateCopy(basis->state(i), gBuffer_.array() + sliceShift);
    sliceShift += stateSize;
  } 
 
  // Perform global communication
  
  setStatus(broadcasting);
  log() << "Allgatherv: " << recv_counts[myCpu]  << " / " << basis->stateCount() * stateSize << " -> " << displacements[numCpus-1] + recv_counts[numCpus-1] << " / " << stateSize * totalStateCount;
  timeCommunicator_->allGatherv(gBuffer_.array() + displacements[myCpu], recv_counts[myCpu], gBuffer_.array(), recv_counts, displacements);
  log() << " complete !\n";
  setStatus(available);
 
  // Update received basis
  setReceivedBroadcastBasis(DynamStateBasisWrapper::New(basis->vectorSize(), totalStateCount, gBuffer_.array()).ptr());
}

void
CommManagerImpl::correctionVectorSizeIs(size_t s) {
  setStatus(receiving);
  rBuffer_.sizeIs(2 * s);
  //RecInfo recInfo;
  log() << "Waiting for cpu " << previousCpuRank().value() << "\n";
  timeCommunicator_->recFrom(previousCpuRank().value(), rBuffer_.array(), rBuffer_.size());
  log() << "Received " << rBuffer_[0] << " from cpu " << previousCpuRank() << "\n";
  setStatus(available);
  setReceivedCorrection(DynamState(s, rBuffer_.array()));
}

void
CommManagerImpl::nextCorrectionIs(const DynamState & state) {
  setStatus(sending);
  sBuffer_.sizeIs(2 * state.vectorSize());
  bufferStateCopy(state, sBuffer_.array());
  log() << "Sending " << sBuffer_[0] << " to cpu " << nextCpuRank() << "\n";
  timeCommunicator_->sendTo(nextCpuRank().value(), localCpuRank().value(), sBuffer_.array(), sBuffer_.size());
  timeCommunicator_->waitForAllReq();
  log() << "Sent " << sBuffer_[0] << " to cpu " << nextCpuRank() << "\n";
  setStatus(available);
  setNextCorrection(state);
}

} // end namespace Pita
