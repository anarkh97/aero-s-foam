#include "SliceMapping.h"

#include <algorithm>

namespace Pita { namespace Hts {

SliceMapping::SliceMapping(FullSliceCount totalFullSlices, CpuCount availableCpus,
                           TaskCount maxWorkload, const SliceStrategy * strategy) :
  taskManager_(TaskManager::New(2 * totalFullSlices.value(), availableCpus.value(), maxWorkload)),
  strategy_(strategy)
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
SliceMapping::hostCpu(SliceId slice) const {
  TaskRank task = strategy_->task(slice);
  WorkerRank worker = taskManager_->worker(task);
  return CpuRank(worker);
}

SliceMapping::SliceIdIterator
SliceMapping::hostedSlice(CpuRank cpu, HalfSliceRank begin, HalfSliceRank end) const {
  SliceIdIterator::ContainerImpl::Ptr containerImpl = new SliceIdIterator::ContainerImpl();
  SliceMapping::SliceIdIterator::SliceIdContainer & slice = containerImpl->slice;

  WorkerRank worker = cpu.value();
  TaskRank beginTask = begin.value();
  TaskRank endTask = end.value();

  for (TaskManager::TaskIterator it = taskManager_->tasks(worker); it; ++it) {
    TaskRank current = *it;
    if (current < beginTask)
      continue;
    if (current >= endTask)
      break;
    
    for (SliceStrategy::SliceIdIterator jt = strategy_->slice(current); jt; ++jt) {
      slice.push_back(*jt);
    }
  }
  
  std::sort(slice.begin(), slice.end()); 

  return SliceIdIterator(containerImpl.ptr());
}


} /* end namespace Hts */ } /* end namespace Pita */
