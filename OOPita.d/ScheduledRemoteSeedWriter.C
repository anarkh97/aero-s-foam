#include "ScheduledRemoteSeedWriter.h"

namespace Pita {
  
ScheduledRemoteSeedWriter::ScheduledRemoteSeedWriter(Communicator * timeComm) :
  RemoteSeedWriter(timeComm),
  phase_(),
  receiveReactor_(NULL)
{}

void
ScheduledRemoteSeedWriter::phaseIs(PhaseRank p) {
  setPhase(p);
  receiveReactor_->notifier()->phaseIs(p);
}

ScheduledRemoteSeedWriter::Manager::Manager(Communicator * timeComm, ScheduledRemoteSeedWriter::Scheduler * scheduler) :
  timeComm_(timeComm),
  scheduler_(scheduler)
{}

void
ScheduledRemoteSeedWriter::ReceiveReactor::onStatus() {
  switch (notifier()->status()) {
    case Activity::executing:
      target()->statusIs(ScheduledRemoteSeedWriter::BUSY);
      break;
    default:
      break;
  }
}

ScheduledRemoteSeedWriter *
ScheduledRemoteSeedWriter::Manager::instanceNew(const Hs::SeedId & key) {
  ScheduledRemoteSeedWriter * i = Impl::instanceNew(key);
  if (scheduler_) {
    scheduler_->instanceNew(key, i);
  }
  return i;
}

void
ScheduledRemoteSeedWriter::Manager::instanceDel(const Hs::SeedId & key) {
  if (scheduler_) {
    scheduler_->instanceDel(key);
  }
  Impl::instanceDel(key);
}


ScheduledRemoteSeedWriter *
ScheduledRemoteSeedWriter::Manager::createNewInstance(const Hs::SeedId & key) {
  return new ScheduledRemoteSeedWriter(timeComm_);
}

} // namespace Pita
