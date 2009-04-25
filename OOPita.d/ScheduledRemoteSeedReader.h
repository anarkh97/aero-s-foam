#ifndef PITA_SCHEDULEDREMOTESEEDREADER_H
#define PITA_SCHEDULEDREMOTESEEDREADER_H

#include "RemoteSeedReader.h"

#include "Activity.h"

namespace Pita {

class ScheduledRemoteSeedReader : public RemoteSeedReader {
public:
  EXPORT_PTRINTERFACE_TYPES(ScheduledRemoteSeedReader);
  
  template <typename Key> class Manager;

  virtual void onState(); // Overriden

  PhaseRank phase() const { return phase_; }
  void phaseIs(PhaseRank p);

protected:
  template <typename Key> friend class Manager;

  class SendReactor : public Activity::Notifiee {
  public:
    typedef Fwk::Ptr<SendReactor> Ptr;
    typedef Fwk::Ptr<const SendReactor> PtrConst;

    virtual void onStatus(); // Overriden

    ScheduledRemoteSeedReader * parent() const { return parent_; }

    SendReactor(Activity * notifier, ScheduledRemoteSeedReader * parent) :
      Activity::Notifiee(notifier),
      parent_(parent)
    {}

  private:
    ScheduledRemoteSeedReader * parent_;
  };

  explicit ScheduledRemoteSeedReader(Communicator * comm);

private:
  PhaseRank phase_;
  SendReactor::Ptr sendReactor_;
};

template <typename Key>
class ScheduledRemoteSeedReader::Manager : public Fwk::PtrInterface<Manager<Key> >, private Fwk::GenManagerImpl<ScheduledRemoteSeedReader, Key> {
public:
  EXPORT_PTRINTERFACE_TYPES(Manager<Key>);

  ScheduledRemoteSeedReader * instance(const Key & key) const { return Fwk::GenManagerImpl<ScheduledRemoteSeedReader, Key>::instance(key); } 
  size_t instanceCount() const { return Fwk::GenManagerImpl<ScheduledRemoteSeedReader, Key>::instanceCount(); }

  ScheduledRemoteSeedReader * instanceNew(const Key & key) { return Fwk::GenManagerImpl<ScheduledRemoteSeedReader, Key>::instanceNew(key); }
  void instanceDel(const Key & key) { Fwk::GenManagerImpl<ScheduledRemoteSeedReader, Key>::instanceDel(key); }

  static Ptr New(Communicator * comm) {
    return new Manager(comm);
  } 

protected:
  explicit Manager(Communicator * comm) :
    timeComm_(comm)
  {}

  virtual ScheduledRemoteSeedReader * createNewInstance(const Key & key);

private:
  Communicator * timeComm_;
};

/* Template fix */
template <typename Key>
ScheduledRemoteSeedReader *
ScheduledRemoteSeedReader::Manager<Key>::createNewInstance(const Key & key) {
  Activity * activity = activityManagerInstance()->activityNew("RemoteReader_" + toString(key)).ptr();
  activity->iterationIs(activityManagerInstance()->currentIteration());
  ScheduledRemoteSeedReader * reader = new ScheduledRemoteSeedReader(timeComm_);
  SendReactor * reactor = new SendReactor(activity, reader);
  reader->sendReactor_ = reactor;
  return reader;
}

} // end namespace Pita

#endif /* PITA_SCHEDULEDREMOTESEEDREADER_H */
