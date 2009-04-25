#ifndef PITA_TIMESLICEMAPPING_H
#define PITA_TIMESLICEMAPPING_H

#include <Utils.d/Connectivity.h>
#include "Fwk.h"
#include "Types.h"
#include <algorithm>

namespace Pita {

class TimeSliceMapping : public PtrInterface<TimeSliceMapping> {
public:
  typedef Fwk::Ptr<TimeSliceMapping> Ptr;
  typedef Fwk::Ptr<const TimeSliceMapping> PtrConst;

  class SliceIteratorConst {
  public:
    SliceIteratorConst & operator++();
    SliceIteratorConst operator++(int) { SliceIteratorConst it(*this); this->operator++(); return it; }
    SliceRank operator*() const { return (mapping_ ? SliceRank(*(mapping_->cpuToTs->operator[](cpuRank_.value()) + offset_.value())) : SliceRank(-1)); }
    operator bool() const { return mapping_; } 
  private:
    friend class TimeSliceMapping;
    SliceIteratorConst(const TimeSliceMapping * mapping, CpuRank cpuRank);
    TimeSliceMapping::PtrConst mapping_;
    CpuRank cpuRank_;
    CpuCount offset_;
  };
  friend class SliceIteratorConst;
  
  // Mutators
  void convergedSlicesInc(SliceCount increment = SliceCount(1));

  // Read Accessors
  SliceCount totalSlices() const { return totalSlices_; }
  SliceCount maxActiveSlices() const { return maxActiveSlices_; }
  CpuCount availableCpus() const { return availableCpus_; }
  
  SliceCount convergedSlices() const { return convergedSlices_; } 
  SliceRank firstActiveSlice() const { return SliceRank(convergedSlices_.value()); }
  SliceRank firstInactiveSlice() const { return firstInactiveSlice_; }
  SliceCount activeSlices() const { return SliceCount(firstInactiveSlice_.value()) - convergedSlices_; }
  
  SliceCount slicesOnCpu(CpuRank cpuRank) const { return SliceCount(cpuToTs->num(cpuRank.value())); } 
  SliceCount activeSlicesOnCpu(CpuRank cpuRank) const;
  CpuRank owningCpu(SliceRank sliceRank) const { return CpuRank(tsToCpu->getTargetValue(sliceRank.value())); } 
  
  SliceIteratorConst slices(CpuRank cpuRank) const { return SliceIteratorConst(this, cpuRank); }

  //bool hasFirstInactiveSlice(int cpuRank) const;

  class Notifiee : public Fwk::BaseNotifiee<TimeSliceMapping> {
  public:
    typedef Fwk::Ptr<Notifiee> Ptr;
    typedef Fwk::Ptr<const Notifiee> PtrConst;

    virtual void onConvergedSlices() {}

  protected:
    explicit Notifiee(TimeSliceMapping * notifier) : Fwk::BaseNotifiee<TimeSliceMapping>(notifier) {}
  };
  
  Notifiee::Ptr lastNotifiee() const { return notifiee_; }
  virtual void lastNotifieeIs(Notifiee * lastNotifiee) { notifiee_ = lastNotifiee; }
  
  static TimeSliceMapping::Ptr New(SliceCount totalSlices, SliceCount maxActiveSlices, CpuCount availableCpus) {
    return new TimeSliceMapping(totalSlices, maxActiveSlices, availableCpus);
  }
  
protected:
  TimeSliceMapping(SliceCount totalSlices, SliceCount maxActiveSlices, CpuCount availableCpus);
  TimeSliceMapping(const TimeSliceMapping &); // No implementation
  TimeSliceMapping & operator=(const TimeSliceMapping &); // No implementation
  ~TimeSliceMapping();
  
private:
  void updateFirstInactiveSlice() { firstInactiveSlice_ = SliceRank(std::min(totalSlices().value(), maxActiveSlices().value() * availableCpus().value() + convergedSlices().value())); }
  
  // Structural variables
  SliceCount totalSlices_;
  SliceCount maxActiveSlices_;
  CpuCount availableCpus_;
 
  // Computational status
  SliceCount convergedSlices_;  // Also corresponds to first active slice number
  SliceRank firstInactiveSlice_;

  // Mapping slices <-> Cpus
  Connectivity * tsToCpu;
  Connectivity * cpuToTs;

  // Notification support
  Notifiee * notifiee_;
};

OStream & operator<<(OStream & out, const TimeSliceMapping & tsm);

} // end namespace Pita

#endif
