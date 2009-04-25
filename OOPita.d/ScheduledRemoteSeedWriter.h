#ifndef PITA_SCHEDULEDREMOTESEEDWRITER_H
#define PITA_SCHEDULEDREMOTESEEDWRITER_H

#include "RemoteSeedWriter.h"

#include "Activity.h"

namespace Pita {

class ScheduledRemoteSeedWriter : public RemoteSeedWriter {
public:
  EXPORT_PTRINTERFACE_TYPES(ScheduledRemoteSeedWriter);

  class Manager;
  class Scheduler;
  
  class ReceiveReactor : public Activity::Notifiee {
  public:
    typedef Fwk::Ptr<ReceiveReactor> Ptr;
    typedef Fwk::Ptr<const ReceiveReactor> PtrConst;

    virtual void onStatus(); // Overriden

    ScheduledRemoteSeedWriter * target() const { return target_; }

    ReceiveReactor(Activity * notifier, ScheduledRemoteSeedWriter * target) :
      Activity::Notifiee(notifier),
      target_(target)
    {}

  private:
    ScheduledRemoteSeedWriter * target_;
  };

  PhaseRank phase() const { return phase_; }
  void phaseIs(PhaseRank p);

protected:
  friend class Manager;
  friend class Scheduler;
  
  ScheduledRemoteSeedWriter(Communicator * timeComm);

  void setOriginCpu(CpuRank oc) { originCpu_ = oc; }
  void setVectorSize(size_t vs) { vectorSize_ = vs; }
  void setStatus(Status s) { status_ = s; }
  void setTarget(Seed * t) { target_ = t; }
  void setPhase(PhaseRank p) { phase_ = p; }

private:
  CpuRank originCpu_;
  Seed::Ptr target_;
  size_t vectorSize_;
  Status status_;
  PhaseRank phase_;
  ReceiveReactor::Ptr receiveReactor_;

  Communicator * timeComm_;
  SimpleBuffer<double> stateBuffer_;
};


class ScheduledRemoteSeedWriter::Scheduler : public Fwk::PtrInterface<ScheduledRemoteSeedWriter::Scheduler> {
public:
  EXPORT_PTRINTERFACE_TYPES(Scheduler);

  friend class ScheduledRemoteSeedWriter::Manager;

protected:
  virtual void instanceNew(const Hs::SeedId & key, ScheduledRemoteSeedWriter *) = 0;
  virtual void instanceDel(const Hs::SeedId & key) = 0;

  ScheduledRemoteSeedWriter::ReceiveReactor * getReceiveReactor(ScheduledRemoteSeedWriter * w) const {
    return w ? w->receiveReactor_.ptr() : NULL;
  }

  void setReceiveReactor(ScheduledRemoteSeedWriter * w, ScheduledRemoteSeedWriter::ReceiveReactor * r) const {
    w->receiveReactor_ = r;
  }
};


class ScheduledRemoteSeedWriter::Manager : public Fwk::PtrInterface<ScheduledRemoteSeedWriter::Manager>, private Fwk::GenManagerImpl<ScheduledRemoteSeedWriter, Hs::SeedId> {
public:
  EXPORT_PTRINTERFACE_TYPES(Manager);
  typedef Fwk::GenManagerImpl<ScheduledRemoteSeedWriter, Hs::SeedId> Impl;

  ScheduledRemoteSeedWriter * instance(const Hs::SeedId & key) const { return Impl::instance(key); }
  InstanceCount instanceCount() const { return Impl::instanceCount(); }
  
  ScheduledRemoteSeedWriter * instanceNew(const Hs::SeedId & key);
  void instanceDel(const Hs::SeedId & key);

  static Ptr New(Communicator * timeComm, ScheduledRemoteSeedWriter::Scheduler * scheduler) { 
    return new Manager(timeComm, scheduler);
  }

protected:
  Manager(Communicator * timeComm, ScheduledRemoteSeedWriter::Scheduler * scheduler);

  virtual ScheduledRemoteSeedWriter * createNewInstance(const Hs::SeedId & key);

private:
  Communicator * timeComm_;
  ScheduledRemoteSeedWriter::Scheduler::Ptr scheduler_;
};
  

} // end namespace Pita

#endif /* PITA_SCHEDULEDREMOTESEEDWRITER_H */
