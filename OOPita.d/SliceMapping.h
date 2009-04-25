#ifndef PITA_HTS_SLICEMAPPING_H
#define PITA_HTS_SLICEMAPPING_H

#include "Fwk.h"
#include "Types.h"

#include "TaskManager.h"
#include "SliceStrategy.h"

#include <vector>

namespace Pita { namespace Hts {

class SliceMapping : public Fwk::PtrInterface<SliceMapping> {
public:
  EXPORT_PTRINTERFACE_TYPES(SliceMapping);

  class SliceIdIterator;

  HalfSliceCount totalSlices() const;
  CpuCount availableCpus() const;
  TaskCount maxWorkload() const;

  HalfSliceCount activeSlices() const;
  HalfSliceRank firstActiveSlice() const;
  HalfSliceRank firstInactiveSlice() const;

  HalfSliceCount convergedSlices() const;
  void convergedSlicesInc(HalfSliceCount increment = HalfSliceCount(1));

  CpuRank hostCpu(SliceId slice) const;
 
  SliceIdIterator hostedSlice(CpuRank cpu,
                              HalfSliceRank begin,
                              HalfSliceRank end) const; 

  static Ptr New(FullSliceCount totalFullSlices, CpuCount availableCpus,
                 TaskCount maxWorkload, const SliceStrategy * strategy) {
    return new SliceMapping(totalFullSlices, availableCpus, maxWorkload, strategy); 
  }

protected:
  SliceMapping(FullSliceCount totalFullSlices, CpuCount availableCpus,
               TaskCount maxWorkload, const SliceStrategy * strategy);

  friend class SliceIdIterator;

private:
  TaskManager::Ptr taskManager_;
  SliceStrategy::PtrConst strategy_;

  DISALLOW_COPY_AND_ASSIGN(SliceMapping);
};


class SliceMapping::SliceIdIterator {
public:
  const SliceId & operator*() const { return *current_; }
  const SliceId * operator->() const { return current_.operator->(); }
  SliceIdIterator & operator++() { ++current_; return *this; }
  SliceIdIterator operator++(int) { SliceIdIterator tmp(*this); ++(*this); return tmp; }
  operator bool() const { return current_ != end_; }

  SliceIdIterator() :
    container_(new ContainerImpl()),
    current_(container_->slice.begin()),
    end_(container_->slice.end())
  {}

  // Default copy, assignment

protected:
  typedef std::vector<SliceId> SliceIdContainer;
  typedef SliceIdContainer::const_iterator ItImpl;
  
  class ContainerImpl : public Fwk::PtrInterface<ContainerImpl> {
  public:
    EXPORT_PTRINTERFACE_TYPES(ContainerImpl);
    ContainerImpl() : slice() {}
    SliceIdContainer slice;
  };

  explicit SliceIdIterator(const ContainerImpl * container) :
    container_(container),
    current_(container_->slice.begin()),
    end_(container_->slice.end())
  {}

  friend class SliceMapping;

private:
  ContainerImpl::PtrConst container_;
  ItImpl current_, end_;  
};

} /* end namespace Hts */
} /* end namespace Pita */

#endif /* PITA_HTS_SLICEMAPPING_H */
