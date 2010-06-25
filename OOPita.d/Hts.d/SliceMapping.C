#include "SliceMapping.h"

#include <algorithm>

namespace Pita { namespace Hts {

SliceMapping::SliceMapping(FullSliceCount totalFullSlices, CpuCount availableCpus, HalfSliceCount maxWorkload) :
  loadBalancer_(LoadBalancer::New(2 * totalFullSlices.value(), availableCpus.value(), maxWorkload.value()))
{}

HalfSliceCount
SliceMapping::totalSlices() const {
  return HalfSliceCount(loadBalancer_->totalTasks());
}

CpuCount
SliceMapping::availableCpus() const {
  return CpuCount(loadBalancer_->availableWorkers());
}


HalfSliceCount
SliceMapping::maxWorkload() const {
  return HalfSliceCount(loadBalancer_->maxWorkload());
}

HalfSliceCount
SliceMapping::activeSlices() const {
  return HalfSliceCount(loadBalancer_->currentGlobalWorkload());
}

HalfSliceRank
SliceMapping::firstActiveSlice() const {
  return HalfSliceRank(loadBalancer_->firstCurrentTask());
}

HalfSliceRank
SliceMapping::firstInactiveSlice() const {
  return HalfSliceRank(loadBalancer_->firstWaitingTask());
}

HalfSliceCount
SliceMapping::convergedSlices() const {
  return HalfSliceCount(loadBalancer_->completedTasks());
}

void
SliceMapping::SliceMapping::convergedSlicesInc(HalfSliceCount increment) {
  loadBalancer_->completedTasksInc(increment.value());
}

CpuRank
SliceMapping::hostCpu(HalfSliceRank slice) const {
  return CpuRank(loadBalancer_->worker(slice.value()));
}

SliceMapping::SliceIterator
SliceMapping::hostedSlice(CpuRank cpu) const {
  return SliceIterator(loadBalancer_->tasks(cpu.value()));
}


} /* end namespace Hts */ } /* end namespace Pita */
