#ifndef PITA_SCHEDULEDREMOTESHAREDSTATE_H
#define PITA_SCHEDULEDREMOTESHAREDSTATE_H

#include "RemoteSharedState.h"

#include "Activity.h"

namespace Pita {

template <typename S>
class ScheduledRemoteSharedStateReader : public RemoteSharedStateReader<S> {
public:
  EXPORT_PTRINTERFACE_TYPES(ScheduledRemoteSharedStateReader);
  
  template <typename Key> class Manager;

  virtual void onState(); // Overriden

  PhaseRank phase() const { return phase_; }
  void phaseIs(PhaseRank p);

protected:
  template <typename Key> friend class Manager;

  class SendReactor : public Activity::Notifiee {
  public:
    EXPORT_PTRINTERFACE_TYPES(SendReactor);

    virtual void onStatus(); // Overriden

    ScheduledRemoteSharedStateReader<S> * parent() const { return parent_; }

    SendReactor(Activity * notifier, ScheduledRemoteSharedStateReader<S> * parent) :
      Activity::Notifiee(notifier),
      parent_(parent)
    {}

  private:
    ScheduledRemoteSharedStateReader<S> * parent_;
  };

  explicit ScheduledRemoteSharedStateReader(Communicator * comm);

private:
  PhaseRank phase_;
  Fwk::Ptr<SendReactor> sendReactor_;
};

template <typename S>
template <typename Key>
class ScheduledRemoteSharedStateReader<S>::Manager : public Fwk::PtrInterface<Manager<Key> >, private Fwk::GenManagerImpl<ScheduledRemoteSharedStateReader<S>, Key> {
public:
  EXPORT_PTRINTERFACE_TYPES(Manager<Key>);

  ScheduledRemoteSharedStateReader<S> * instance(const Key & key) const { return Fwk::GenManagerImpl<ScheduledRemoteSharedStateReader<S>, Key>::instance(key); } 
  size_t instanceCount() const { return Fwk::GenManagerImpl<ScheduledRemoteSharedStateReader<S>, Key>::instanceCount(); }

  ScheduledRemoteSharedStateReader<S> * instanceNew(const Key & key) { return Fwk::GenManagerImpl<ScheduledRemoteSharedStateReader<S>, Key>::instanceNew(key); }
  void instanceDel(const Key & key) { Fwk::GenManagerImpl<ScheduledRemoteSharedStateReader<S>, Key>::instanceDel(key); }

  static Ptr New(Communicator * comm) {
    return new Manager(comm);
  } 

protected:
  explicit Manager(Communicator * comm) :
    timeComm_(comm)
  {}

  virtual ScheduledRemoteSharedStateReader<S> * createNewInstance(const Key & key);

private:
  Communicator * timeComm_;
};

template <typename S>
template <typename Key>
ScheduledRemoteSharedStateReader<S> *
ScheduledRemoteSharedStateReader<S>::Manager<Key>::createNewInstance(const Key & key) {
  Activity * activity = activityManagerInstance()->activityNew("RemoteReader_" + toString(key)).ptr();
  activity->iterationIs(activityManagerInstance()->currentIteration());
  ScheduledRemoteSharedStateReader<S> * reader = new ScheduledRemoteSharedStateReader<S>(timeComm_);
  SendReactor * reactor = new SendReactor(activity, reader);
  reader->sendReactor_ = reactor;
  return reader;
}


template <typename S>
ScheduledRemoteSharedStateReader<S>::ScheduledRemoteSharedStateReader(Communicator * comm) :
  RemoteSharedStateReader<S>(comm),
  phase_(),
  sendReactor_(NULL)
{}

template <typename S>
void
ScheduledRemoteSharedStateReader<S>::onState() {
  if (this->notifier()->status() != SharedState<S>::INACTIVE && this->status() == ScheduledRemoteSharedStateReader<S>::READY) {
    statusIs(ScheduledRemoteSharedStateReader<S>::SCHEDULED);
    Activity * a = sendReactor_->notifier().ptr();
    IterationRank iteration = a->currentPhase() <= a->phase() ? a->currentIteration() : a->currentIteration().next();
    a->iterationIs(iteration);
    a->statusIs(Activity::scheduled);
  }
}

template <typename S>
void
ScheduledRemoteSharedStateReader<S>::phaseIs(PhaseRank p) {
  phase_ = p;
  sendReactor_->notifier()->phaseIs(p);
}

template <typename S>
void
ScheduledRemoteSharedStateReader<S>::SendReactor::onStatus() {
  switch (notifier()->status()) {
    case Activity::executing:
      parent()->statusIs(ScheduledRemoteSharedStateReader<S>::BUSY);
      break;
    default:
      break;
  }
}

/* Writer */

template <typename S>
class ScheduledRemoteSharedStateWriter : public RemoteSharedStateWriter<S> {
public:
  EXPORT_PTRINTERFACE_TYPES(ScheduledRemoteSharedStateWriter);

  class Manager;
  class Scheduler;
  
  class ReceiveReactor : public Activity::Notifiee {
  public:
    typedef Fwk::Ptr<ReceiveReactor> Ptr;
    typedef Fwk::Ptr<const ReceiveReactor> PtrConst;

    virtual void onStatus(); // Overriden

    ScheduledRemoteSharedStateWriter<S> * target() const { return target_; }

    ReceiveReactor(Activity * notifier, ScheduledRemoteSharedStateWriter<S> * target) :
      Activity::Notifiee(notifier),
      target_(target)
    {}

  private:
    ScheduledRemoteSharedStateWriter<S> * target_;
  };

  PhaseRank phase() const { return phase_; }
  void phaseIs(PhaseRank p);

protected:
  friend class Manager;
  friend class Scheduler;
  
  ScheduledRemoteSharedStateWriter(Communicator * timeComm);

  void setOriginCpu(CpuRank oc) { originCpu_ = oc; }
  void setVectorSize(size_t vs) { vectorSize_ = vs; }
  void setStatus(typename RemoteSharedStateWriter<S>::Status s) { status_ = s; }
  void setTarget(SharedState<S> * t) { target_ = t; }
  void setPhase(PhaseRank p) { phase_ = p; }

private:
  CpuRank originCpu_;
  Fwk::Ptr<SharedState<S> > target_;
  size_t vectorSize_;
  typename RemoteSharedStateWriter<S>::Status status_;
  PhaseRank phase_;
  Fwk::Ptr<ReceiveReactor> receiveReactor_;

  Communicator * timeComm_;
  SimpleBuffer<double> stateBuffer_;
};

template <typename S>
class ScheduledRemoteSharedStateWriter<S>::Scheduler : public Fwk::PtrInterface<Scheduler> {
public:
  EXPORT_PTRINTERFACE_TYPES(Scheduler);

  friend class ScheduledRemoteSharedStateWriter::Manager;

protected:
  virtual void instanceNew(const Hs::SeedId & key, ScheduledRemoteSharedStateWriter<S> *) = 0;
  virtual void instanceDel(const Hs::SeedId & key) = 0;

  typename ScheduledRemoteSharedStateWriter<S>::ReceiveReactor * getReceiveReactor(ScheduledRemoteSharedStateWriter<S> * w) const {
    return w ? w->receiveReactor_.ptr() : NULL;
  }

  void setReceiveReactor(ScheduledRemoteSharedStateWriter<S> * w, typename ScheduledRemoteSharedStateWriter<S>::ReceiveReactor * r) const {
    w->receiveReactor_ = r;
  }
};

template <typename S>
class ScheduledRemoteSharedStateWriter<S>::Manager : public Fwk::PtrInterface<Manager>, private Fwk::GenManagerImpl<ScheduledRemoteSharedStateWriter<S>, Hs::SeedId> {
public:
  EXPORT_PTRINTERFACE_TYPES(Manager);
  typedef Fwk::GenManagerImpl<ScheduledRemoteSharedStateWriter, Hs::SeedId> Impl;

  ScheduledRemoteSharedStateWriter * instance(const Hs::SeedId & key) const { return Impl::instance(key); }
  typename Fwk::GenManagerImpl<ScheduledRemoteSharedStateWriter<S>, Hs::SeedId>::InstanceCount instanceCount() const { return Impl::instanceCount(); }
  
  ScheduledRemoteSharedStateWriter * instanceNew(const Hs::SeedId & key);
  void instanceDel(const Hs::SeedId & key);

  static Ptr New(Communicator * timeComm, ScheduledRemoteSharedStateWriter::Scheduler * scheduler) { 
    return new Manager(timeComm, scheduler);
  }

protected:
  Manager(Communicator * timeComm, ScheduledRemoteSharedStateWriter::Scheduler * scheduler);

  virtual ScheduledRemoteSharedStateWriter * createNewInstance(const Hs::SeedId & key);

private:
  Communicator * timeComm_;
  Fwk::Ptr<typename ScheduledRemoteSharedStateWriter<S>::Scheduler> scheduler_;
};

template <typename S>
ScheduledRemoteSharedStateWriter<S>::ScheduledRemoteSharedStateWriter(Communicator * timeComm) :
  RemoteSharedStateWriter<S>(timeComm),
  phase_(),
  receiveReactor_(NULL)
{}

template <typename S>
void
ScheduledRemoteSharedStateWriter<S>::phaseIs(PhaseRank p) {
  setPhase(p);
  receiveReactor_->notifier()->phaseIs(p);
}

template <typename S>
ScheduledRemoteSharedStateWriter<S>::Manager::Manager(Communicator * timeComm, ScheduledRemoteSharedStateWriter::Scheduler * scheduler) :
  timeComm_(timeComm),
  scheduler_(scheduler)
{}

template <typename S>
void
ScheduledRemoteSharedStateWriter<S>::ReceiveReactor::onStatus() {
  switch (notifier()->status()) {
    case Activity::executing:
      target()->statusIs(ScheduledRemoteSharedStateWriter::BUSY);
      break;
    default:
      break;
  }
}

template <typename S>
ScheduledRemoteSharedStateWriter<S> *
ScheduledRemoteSharedStateWriter<S>::Manager::instanceNew(const Hs::SeedId & key) {
  ScheduledRemoteSharedStateWriter * i = Impl::instanceNew(key);
  if (scheduler_) {
    scheduler_->instanceNew(key, i);
  }
  return i;
}

template <typename S>
void
ScheduledRemoteSharedStateWriter<S>::Manager::instanceDel(const Hs::SeedId & key) {
  if (scheduler_) {
    scheduler_->instanceDel(key);
  }
  Impl::instanceDel(key);
}


template <typename S>
ScheduledRemoteSharedStateWriter<S> *
ScheduledRemoteSharedStateWriter<S>::Manager::createNewInstance(const Hs::SeedId & key) {
  return new ScheduledRemoteSharedStateWriter<S>(timeComm_);
}

} // end namespace Pita

#endif /* PITA_SCHEDULEDREMOTESHAREDSTATE_H */
