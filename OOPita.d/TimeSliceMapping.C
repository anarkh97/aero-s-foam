#include "TimeSliceMapping.h"
#include <algorithm>

namespace Pita {

TimeSliceMapping::TimeSliceMapping(SliceCount nTs, SliceCount nMaxActive, CpuCount nCpu) :
  totalSlices_(nTs),
  maxActiveSlices_(nMaxActive),
  availableCpus_(nCpu),
  convergedSlices_(0),
  firstInactiveSlice_(0),
  tsToCpu(NULL),
  cpuToTs(NULL)
{
  updateFirstInactiveSlice();

  int totalSliceCount = totalSlices().value();
  int availableCpuCount = availableCpus().value();
  int maxActiveSliceCount = maxActiveSlices().value();
  
  int * arrayTs = new int[totalSliceCount + 1];
  for (int i = 0; i <= totalSliceCount; ++i)
    arrayTs[i] = i;

  int * targetCpu = new int[totalSliceCount];
  
  int ratio, remain;
  if (totalSliceCount < availableCpuCount * maxActiveSliceCount) {
    ratio = totalSliceCount / availableCpuCount;
    remain = totalSliceCount - ratio * availableCpuCount;
  } else {
    ratio = maxActiveSliceCount;
    remain = 0;
  }

  int currentCpu = 0;
  int remainingSpotsOnCurrentCpu = ratio + (remain > 0 ? 1 : 0);
  for (int currentTs = 0; currentTs < totalSliceCount; ++currentTs) {
    if (remainingSpotsOnCurrentCpu <= 0) {
      currentCpu = (currentCpu + 1) % availableCpuCount;
      remainingSpotsOnCurrentCpu = ratio + (remain > currentCpu ? 1 : 0);
    }
    targetCpu[currentTs] = currentCpu;
    --remainingSpotsOnCurrentCpu;
  }

  tsToCpu = new Connectivity(totalSliceCount, arrayTs, targetCpu, 1);
  cpuToTs = tsToCpu->reverse();
  cpuToTs->sortTargets(); 
}

TimeSliceMapping::~TimeSliceMapping() {
  delete cpuToTs;
  delete tsToCpu;
}

void
TimeSliceMapping::convergedSlicesInc(SliceCount increment) {
  if (increment != SliceCount(0)) {
    convergedSlices_ += increment;
    updateFirstInactiveSlice();
    if (notifiee_) {
      notifiee_->onConvergedSlices();
    }
  }
}

SliceCount
TimeSliceMapping::activeSlicesOnCpu(CpuRank cpuRank) const {
  SliceCount activeSliceCount(0);
  for (SliceIteratorConst it = slices(cpuRank); it && (*it < firstInactiveSlice()); ++it) {
    if ((*it).value() >= convergedSlices().value())
      ++activeSliceCount;
  }
  return activeSliceCount;
}

/*bool TimeSliceMapping::hasFirstInactiveSlice(int Cpuid) const
{
  if (firstInactiveSlice_ >= numSlices_)
    return false;
  return (*(tsToCpu->operator[](firstInactiveSlice_)) == Cpuid);
}*/

TimeSliceMapping::SliceIteratorConst::SliceIteratorConst(const TimeSliceMapping * mapping, CpuRank cpuRank) : mapping_(mapping), cpuRank_(cpuRank), offset_(0u) {
  if (cpuRank.value() >= std::min(mapping_->totalSlices().value(), mapping_->availableCpus().value()))
    mapping_ = NULL;    
}

TimeSliceMapping::SliceIteratorConst &
TimeSliceMapping::SliceIteratorConst::operator++() {
  if (mapping_) {
    if ((++offset_).value() >= mapping_->slicesOnCpu(cpuRank_).value())
       mapping_ = NULL;
  }
  return *this;
}

OStream & operator<<(OStream & out, const TimeSliceMapping & tsm) {
  for (CpuRank cpu(0); cpu.value() < tsm.availableCpus().value(); ++cpu) {
    out << "Cpu # " << cpu << " ->";
    for (TimeSliceMapping::SliceIteratorConst it = tsm.slices(cpu); it; ++it) {
      out << ' ' << *it;
    }
    out << '\n';
  }
 return out; 
}

} // end namespace Pita
