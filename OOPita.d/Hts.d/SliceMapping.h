#ifndef PITA_HTS_SLICEMAPPING_H
#define PITA_HTS_SLICEMAPPING_H

#include "Fwk.h"
#include "Types.h"

#include "../LoadBalancer.h"

#include <vector>

namespace Pita { namespace Hts {

class SliceMapping : public Fwk::PtrInterface<SliceMapping> {
public:
  EXPORT_PTRINTERFACE_TYPES(SliceMapping);

  class SliceIterator;

  HalfSliceCount totalSlices() const;
  CpuCount availableCpus() const;
  TaskCount maxWorkload() const;

  HalfSliceCount activeSlices() const;
  HalfSliceRank firstActiveSlice() const;
  HalfSliceRank firstInactiveSlice() const;

  HalfSliceCount convergedSlices() const;
  void convergedSlicesInc(HalfSliceCount increment = HalfSliceCount(1));

  CpuRank hostCpu(HalfSliceRank slice) const;
 
  SliceIterator hostedSlice(CpuRank cpu) const;

  static Ptr New(FullSliceCount totalFullSlices, CpuCount availableCpus, TaskCount maxWorkload) {
    return new SliceMapping(totalFullSlices, availableCpus, maxWorkload); 
  }

protected:
  SliceMapping(FullSliceCount totalFullSlices, CpuCount availableCpus, TaskCount maxWorkload);

  friend class SliceIterator;

private:
  LoadBalancer::Ptr taskManager_;

  DISALLOW_COPY_AND_ASSIGN(SliceMapping);
};


class SliceMapping::SliceIterator {
public:
  HalfSliceRank operator*() const { return HalfSliceRank(*impl_); }
  SliceIterator & operator++() { ++impl_; return *this; }
  SliceIterator operator++(int) { SliceIterator tmp(*this); ++(*this); return tmp; }
  operator bool() const { return impl_; }

  explicit SliceIterator(LoadBalancer::TaskIterator impl) :
    impl_(impl)
  {}

  // Default copy, assignment

private:
  LoadBalancer::TaskIterator impl_;
};

} /* end namespace Hts */ } /* end namespace Pita */

#endif /* PITA_HTS_SLICEMAPPING_H */
