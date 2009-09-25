#ifndef PITA_HTS_SLICESTRATEGY_H
#define PITA_HTS_SLICESTRATEGY_H

#include "Fwk.h"
#include "Types.h"

#include "../TaskManager.h"

#include <vector>

namespace Pita { namespace Hts {

class SliceStrategy : public Fwk::PtrInterface<SliceStrategy> {
public:
  EXPORT_PTRINTERFACE_TYPES(SliceStrategy);

  class SliceIdIterator;

  virtual TaskRank task(const SliceId & id) const = 0;
  virtual SliceIdIterator slice(TaskRank workLoadId) const = 0;

protected:
  friend class SliceIdIterator;

  SliceStrategy() {}

  SliceIdIterator iteratorNew(const std::vector<SliceId> & slice) const;

  DISALLOW_COPY_AND_ASSIGN(SliceStrategy);
};


class SliceStrategy::SliceIdIterator {
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

  friend class SliceStrategy;

private:
  ContainerImpl::PtrConst container_;
  ItImpl current_, end_;
};

inline
SliceStrategy::SliceIdIterator
SliceStrategy::iteratorNew(const SliceStrategy::SliceIdIterator::SliceIdContainer & slice) const {
  SliceIdIterator::ContainerImpl::Ptr container = new SliceIdIterator::ContainerImpl();
  container->slice = slice;
  return SliceIdIterator(container.ptr());
}

} // end namespace Hts
} // end namespace Pita

#endif /* PITA_HTS_SLICESTRATEGY_H */
