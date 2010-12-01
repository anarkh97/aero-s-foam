#include "GlobalStateSharing.h"

#include "../DynamStatePlainBasis.h"
#include "../DynamStateBasisWrapper.h"
#include <Comm.d/Communicator.h>

#include <numeric>

namespace Pita { namespace Hts {

GlobalStateSharing::GlobalStateSharing(Communicator * timeComm, size_t vectorSize) :
  NamedTask("Exchange Global Projection Basis"),
  timeComm_(timeComm),
  vectorSize_(vectorSize),
  seedGetter_(NULL),
  stateCount_(0),
  localStates_(),
  buffer_(),
  bufferCounts_(),
  bufferStrides_(),
  consolidatedBasis_(DynamStatePlainBasis::New(vectorSize))
{}

inline
CpuRank
GlobalStateSharing::localCpu() const {
  return CpuRank(timeComm_->myID());
}

void
GlobalStateSharing::mappingIs(const SliceMapping & mapping) {
  // Allocate buffer
  stateCount_ = mapping.activeSlices().value();

  const int stateSize = 2 * vectorSize();
  const size_t newBufferSize = stateCount_ * stateSize;
  if (buffer_.size() < newBufferSize) {
    buffer_.sizeIs(newBufferSize);
  }
 
  // Partition buffer
  const int cpuCount = mapping.availableCpus().value();
  bufferCounts_.sizeIs(cpuCount);
  bufferStrides_.sizeIs(cpuCount);

  for (int cpu = 0; cpu < cpuCount; ++cpu) {
    bufferCounts_[cpu] = mapping.activeSlices(CpuRank(cpu)).value() * stateSize;
  }

  bufferStrides_[0] = 0;
  std::partial_sum(bufferCounts_.array(),
                   bufferCounts_.array() + cpuCount - 1,
                   bufferStrides_.array() + 1);

  // Enqueue local seeds to be exchanged
  assert(localStates_.empty());
  
  for (SliceMapping::SliceIterator it = mapping.hostedSlice(localCpu()); it; ++it) {
    HalfSliceRank rank = *it;

    if (rank < mapping.firstActiveSlice()) {
      continue;
    }
    if (rank >= mapping.firstInactiveSlice()) {
      break;
    }

    HalfSliceCount distance = rank - mapping.firstActiveSlice();

    SeedId id = (distance.value() % 2) ? SeedId(LEFT_SEED, rank.next()) : SeedId(RIGHT_SEED, rank);
    localStates_.push(seedGetter_(id));
  }
}

void
GlobalStateSharing::iterationIs(IterationRank iter) {
  // Fill local buffer
  const int myCpu = localCpu().value();
  double * bufferBegin = buffer_.array() + bufferStrides_[myCpu];
  const int stateSize = 2 * vectorSize();

  assert(localStates_.size() * stateSize == bufferCounts_[myCpu]);

  while (!localStates_.empty()) {
    const DynamState s = localStates_.front()->state();
    bufferStateCopy(s, bufferBegin);
    bufferBegin += stateSize;
    localStates_.pop();
  }

  // Exchange data
  timeComm_->allGatherv(buffer_.array() + bufferStrides_[myCpu],
                        bufferCounts_[myCpu],
                        buffer_.array(),
                        bufferCounts_.array(),
                        bufferStrides_.array());

  // Expose data
  consolidatedBasis_ = DynamStateBasisWrapper::New(vectorSize(), stateCount_, buffer_.array());
  setIteration(iter);
}

} /* end namespace Hts */ } /* end namespace Pita */
