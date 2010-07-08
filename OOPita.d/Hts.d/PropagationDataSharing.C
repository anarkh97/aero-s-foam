#include "PropagationDataSharing.h"

#include "../DynamStatePlainBasis.h"
#include "../DynamStateBasisWrapper.h"
#include <Comm.d/Communicator.h>

#include <numeric>

namespace Pita { namespace Hts {

void
PropagationDataSharing::iterationIs(IterationRank iter) {
  // Allocate buffer
  const int stateCount = mapping_->activeSlices().value();

  const int stateSize = 2 * vectorSize();
  const size_t newBufferSize = stateCount * stateSize;
  if (buffer_.size() < newBufferSize) {
    buffer_.sizeIs(newBufferSize);
  }
 
  // Partition buffer
  const int cpuCount = mapping_->availableCpus().value();
  SimpleBuffer<int> recvs_counts(cpuCount);
  SimpleBuffer<int> displacements(cpuCount); 

  for (int cpu = 0; cpu < cpuCount; ++cpu) {
    recvs_counts[cpu] = mapping_->activeSlices(CpuRank(cpu)).value() * stateSize;
  }

  displacements[0] = 0;
  std::partial_sum(recvs_counts.array(),
                   recvs_counts.array() + cpuCount - 1,
                   displacements.array() + 1);

  // Fill local buffer
  const int myCpu = timeComm_->myID();
  double * bufferBegin = buffer_.array() + displacements[myCpu];

  // 1) Final states
  BasisCollectorImpl::CollectedState cs = collector_->firstForwardFinalState();
  while (cs.second.vectorSize() != 0) {
    bufferStateCopy(cs.second, bufferBegin);
    bufferBegin += stateSize; 
    collector_->firstForwardFinalStateDel();
    cs = collector_->firstForwardFinalState();
  }

  // 2) Initial states
  cs = collector_->firstBackwardFinalState();
  while (cs.second.vectorSize() != 0) {
    bufferStateCopy(cs.second, bufferBegin);
    bufferBegin += stateSize;
    collector_->firstBackwardFinalStateDel();
    cs = collector_->firstBackwardFinalState();
  }

  // Exchange data
  timeComm_->allGatherv(buffer_.array() + displacements[myCpu],
                        recvs_counts[myCpu],
                        buffer_.array(),
                        recvs_counts.array(),
                        displacements.array());
  // Expose data
  consolidatedBasis_ = DynamStateBasisWrapper::New(vectorSize(), stateCount, buffer_.array());

  setIteration(iter);
}

PropagationDataSharing::PropagationDataSharing(const SliceMapping * mapping, CpuRank localCpu, Communicator * timeComm, size_t vectorSize) :
  NamedTask("Exchange Global Projection Basis"),
  mapping_(mapping),
  localCpu_(localCpu),
  timeComm_(timeComm),
  vectorSize_(vectorSize),
  collector_(BasisCollectorImpl::New()),
  consolidatedBasis_(DynamStatePlainBasis::New(vectorSize))
{}

} /* end namespace Hts */ } /* end namespace Pita */
