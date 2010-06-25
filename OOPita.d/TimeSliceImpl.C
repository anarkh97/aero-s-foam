#include "TimeSliceImpl.h"

namespace Pita {

class TimeSliceImpl::InitialInterface : public TimeSlice::InitialInterface {
public:
  typedef Fwk::Ptr<TimeSliceImpl::InitialInterface> Ptr;
  typedef Fwk::Ptr<const TimeSliceImpl::InitialInterface> PtrConst;

  virtual void seedIs(const DynamState & s) {
    setSeed(s);
    parent_->onSeed();
  }
  virtual void statusIs(Status status) { setStatus(status); }
  
  virtual void onNextSeed() { seedIs(notifier()->nextSeed()); }
  virtual void onStatus() { statusIs(notifier()->status()); }
  
  static InitialInterface::Ptr New(TimeSliceImpl * parent) { return new TimeSliceImpl::InitialInterface(parent); }

protected:
  explicit InitialInterface(TimeSliceImpl * parent) : parent_(parent) {}

private:
  TimeSliceImpl * parent_;
};


class TimeSliceImpl::FinalInterface : public TimeSlice::FinalInterface {
public:
  typedef Fwk::Ptr<TimeSliceImpl::FinalInterface> Ptr;
  typedef Fwk::Ptr<const TimeSliceImpl::FinalInterface> PtrConst;

  static FinalInterface::Ptr New() { return new TimeSliceImpl::FinalInterface(); }

protected:
  FinalInterface() {}
};


// TimesliceImpl Implementation

TimeSliceImpl::TimeSliceImpl(SliceRank rank, Status status) : 
  TimeSlice(rank, status)
{
  setInitialInterface(TimeSliceImpl::InitialInterface::New(this));
  setFinalInterface(TimeSliceImpl::FinalInterface::New());
}

} // end namespace Pita

