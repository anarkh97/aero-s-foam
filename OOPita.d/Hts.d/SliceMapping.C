#include "SliceMapping.h"

#include <algorithm>

namespace Pita { namespace Hts {

SliceMapping::SliceMapping(FullSliceCount totalFullSlices, CpuCount availableCpus, TaskCount maxWorkload) :
  taskManager_(LoadBalancer::New(2 * totalFullSlices.value(), availableCpus.value(), maxWorkload))
{}

HalfSliceCount
SliceMapping::totalSlices() const {
  return HalfSliceCount(taskManager_->totalTasks());
}

CpuCount
SliceMapping::availableCpus() const {
  return CpuCount(taskManager_->availableWorkers());
}


TaskCount
SliceMapping::maxWorkload() const {
  return taskManager_->maxWorkload();
}

HalfSliceCount
SliceMapping::activeSlices() const {
  return HalfSliceCount(taskManager_->currentGlobalWorkload());
}

HalfSliceRank
SliceMapping::firstActiveSlice() const {
  return HalfSliceRank(taskManager_->firstCurrentTask());
}

HalfSliceRank
SliceMapping::firstInactiveSlice() const {
  return HalfSliceRank(taskManager_->firstWaitingTask());
}

HalfSliceCount
SliceMapping::convergedSlices() const {
  return HalfSliceCount(taskManager_->completedTasks());
}

void
SliceMapping::SliceMapping::convergedSlicesInc(HalfSliceCount increment) {
  taskManager_->completedTasksInc(increment.value());
}

CpuRank
SliceMapping::hostCpu(HalfSliceRank slice) const {
  return CpuRank(taskManager_->worker(slice.value()));
}

SliceMapping::SliceIterator
SliceMapping::hostedSlice(CpuRank cpu) const {
  return SliceIterator(taskManager_->tasks(cpu.value()));
}


} /* end namespace Hts */ } /* end namespace Pita */
