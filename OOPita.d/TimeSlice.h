#ifndef PITA_TIMESLICE_H
#define PITA_TIMESLICE_H

#include "Fwk.h"
#include "Types.h"
#include "DynamState.h"
#include "Activity.h"

namespace Pita {

class TimeSlice : public Fwk::PtrInterface<TimeSlice> {
public:
  typedef Fwk::Ptr<TimeSlice> Ptr;
  typedef Fwk::Ptr<const TimeSlice> PtrConst;
  
  SliceRank rank() const { return rank_; };
  
  enum Status {
    inactive = 0,
    active,
    lastIteration,
    converged,
    stopped,
  };
  
  Status status() const { return status_; }
  virtual void statusIs(Status s) = 0;

  class InitialInterface;
  class FinalInterface;
  
  Fwk::Ptr<InitialInterface> initialInterface() { return initialInterface_; }
  Fwk::Ptr<const InitialInterface> initialInterface() const { return initialInterface_; }
  
  Fwk::Ptr<FinalInterface> finalInterface() { return finalInterface_; }
  Fwk::Ptr<const FinalInterface> finalInterface() const { return finalInterface_; }
  
  class Notifiee : public Fwk::BaseNotifiee<TimeSlice> {
  public:
    typedef Fwk::Ptr<Notifiee> Ptr;
    typedef Fwk::Ptr<const Notifiee> PtrConst;

    virtual void onStatus() {}

  protected:
    explicit Notifiee(TimeSlice * notifier) : Fwk::BaseNotifiee<TimeSlice>(notifier) {}
  };
  
  Notifiee::Ptr lastNotifiee() const { return notifiee_; }
  virtual void lastNotifieeIs(Notifiee * notifiee) { setLastNotifiee(notifiee); }
  
  class Manager;
  
  class CorrectionReactor;
  virtual CorrectionReactor * correctionReactor() const = 0;

protected:
  explicit TimeSlice(SliceRank rank, Status status = stopped);
  
  void setStatus(Status newStatus) { status_ = newStatus; if (lastNotifiee()) lastNotifiee()->onStatus(); }
  void setLastNotifiee(Notifiee * notifiee) { notifiee_ = notifiee; }
  void setInitialInterface(Fwk::Ptr<InitialInterface> i) { initialInterface_ = i; }
  void setFinalInterface(Fwk::Ptr<FinalInterface> i) { finalInterface_ = i; }
  
private:
  SliceRank rank_;
  Status status_;
  Fwk::Ptr<InitialInterface> initialInterface_;
  Fwk::Ptr<FinalInterface> finalInterface_;
  Notifiee * notifiee_;
  
  TimeSlice(const TimeSlice &); // No implementation
  TimeSlice & operator=(const TimeSlice &); // No implementation
};


class TimeSlice::FinalInterface : public Fwk::PtrInterface<TimeSlice::FinalInterface> {
public:
  typedef Fwk::Ptr<FinalInterface> Ptr;
  typedef Fwk::Ptr<const FinalInterface> PtrConst;
  
  DynamState nextSeed() const { return nextSeed_; }
  virtual void nextSeedIs(const DynamState & nextSeed) { setNextSeed(nextSeed); }

  enum Status {
    improving = 0,
    final
  };
  
  Status status() const { return status_; }
  void statusIs(Status s) { setStatus(s); }
  
  class Notifiee : public Fwk::BaseNotifiee<FinalInterface> {
  public:
    typedef Fwk::Ptr<Notifiee> Ptr;
    typedef Fwk::Ptr<const Notifiee> PtrConst;

    virtual void onNextSeed() { }
    virtual void onStatus() { }

    explicit Notifiee(FinalInterface * notifier) : BaseNotifiee<FinalInterface>(notifier) {}
  };
  
  Notifiee::Ptr lastNotifiee() const { return notifiee_; } 
  virtual void lastNotifieeIs(const Notifiee::Ptr & notifiee) { setLastNotifiee(notifiee.ptr()); }

protected:
  void setNextSeed(const DynamState & newSeed) { nextSeed_ = newSeed; if (lastNotifiee()) lastNotifiee()->onNextSeed(); }
  void setStatus(Status status) { status_ = status; if (lastNotifiee()) lastNotifiee()->onStatus(); }
  void setLastNotifiee(Notifiee * notifiee) { notifiee_ = notifiee; }
  
  FinalInterface() :
    nextSeed_(),
    status_(improving),
    notifiee_(NULL)
  {}

  FinalInterface(const FinalInterface &); // No implementation
  FinalInterface & operator=(const FinalInterface &); // No implementation
  
private:
  DynamState nextSeed_;
  Status status_;
  Notifiee * notifiee_;
};


class TimeSlice::InitialInterface : public TimeSlice::FinalInterface::Notifiee {
public:
  typedef Fwk::Ptr<InitialInterface> Ptr;
  typedef Fwk::Ptr<const InitialInterface> PtrConst;

  typedef TimeSlice::FinalInterface::Status Status;
  
  DynamState seed() const { return seed_; }
  virtual void seedIs(const DynamState & seed) = 0;

  Status status() const { return status_; }
  virtual void statusIs(Status status) = 0;
  

protected:
  void setSeed(const DynamState & newSeed) { seed_ = newSeed; }  
  void setStatus(Status status) { status_ = status; }
  
  explicit InitialInterface(FinalInterface * notifier = NULL) : FinalInterface::Notifiee(notifier) {}

  InitialInterface(const InitialInterface &); // No implementation
  InitialInterface & operator=(const InitialInterface &); // No implementation
  
private:
  DynamState seed_;
  Status status_;
};


template <typename SliceType>
class TimeSliceIteratorConst {
public:
  explicit TimeSliceIteratorConst(const typename SliceType::Manager * m) :
    sliceManager_(m),
    currentSlice_(m->getNextTimeSlice(NULL))
  {}

  template <typename OtherSliceType>
  TimeSliceIteratorConst(const TimeSliceIteratorConst<OtherSliceType> & it) :
    sliceManager_(it.parent()),
    currentSlice_((*it).ptr())
  {}
  
  template <typename OtherSliceType>
  TimeSliceIteratorConst & operator=(const TimeSliceIteratorConst<OtherSliceType> & it) {
    sliceManager_ = it.parent();
    currentSlice_ = (*it).ptr();
    return *this;
  }

  typename SliceType::Ptr operator*() const { return currentSlice_; } 
  typename SliceType::Ptr operator->() const { return currentSlice_; }

  TimeSliceIteratorConst<SliceType> & operator++();
  TimeSliceIteratorConst<SliceType> operator++(int);
  
  operator bool() const { return currentSlice_ != NULL; }

  const typename SliceType::Manager * parent() const { return sliceManager_; }

private:
  const typename SliceType::Manager * sliceManager_;
  SliceType * currentSlice_;
};

template <typename SliceType>
TimeSliceIteratorConst<SliceType> &
TimeSliceIteratorConst<SliceType>::operator++() {
  if (currentSlice_)
    currentSlice_ = sliceManager_->getNextTimeSlice(currentSlice_);
  return *this;
}

template <typename SliceType>
TimeSliceIteratorConst<SliceType>
TimeSliceIteratorConst<SliceType>::operator++(int) {
  TimeSliceIteratorConst<SliceType> temp(*this);
  ++(*this);
  return temp;
}


class TimeSlice::Manager : public Fwk::PtrInterface<TimeSlice::Manager> {
public:
  typedef Fwk::Ptr<Manager> Ptr;
  typedef Fwk::Ptr<const Manager> PtrConst;

  TimeSlice::Ptr timeSlice(SliceRank rank) const { return findTimeSlice(rank); }
  virtual SliceCount timeSliceCount() const = 0;
  
  TimeSlice::Ptr timeSliceNew(SliceRank rank) { return createTimeSlice(rank); }
  virtual void timeSliceDel(SliceRank rank) = 0;

  friend class TimeSliceIteratorConst<TimeSlice>;
  TimeSliceIteratorConst<TimeSlice> timeSliceIteratorConst() const {
    return TimeSliceIteratorConst<TimeSlice>(this);
  }
  
protected:
  Manager() {}

  // Implementation functions to enable covariant return type
  virtual TimeSlice * createTimeSlice(SliceRank rank) = 0;
  virtual TimeSlice * findTimeSlice(SliceRank rank) const = 0;
  
  // Implementation function for iterator
  virtual TimeSlice * getNextTimeSlice(const TimeSlice *) const = 0;
  
private:
  Manager(Manager &); // No implementation
  Manager & operator=(const Manager &); // No implementation
};

class TimeSlice::CorrectionReactor : public Activity::Notifiee {
public:
  typedef Fwk::Ptr<CorrectionReactor> Ptr;
  typedef Fwk::Ptr<const CorrectionReactor> PtrConst;
  
  virtual const TimeSlice * slice() const = 0;

  explicit CorrectionReactor(Activity * notifier) :
    Activity::Notifiee(notifier)
  {}
};


} // end namespace Pita

#endif /* PITA_TIMESLICE_H */
