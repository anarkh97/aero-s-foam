#include "ScheduledRemoteSeedReader.h"

namespace Pita {

ScheduledRemoteSeedReader::ScheduledRemoteSeedReader(Communicator * comm) :
  RemoteSeedReader(comm),
  phase_(),
  sendReactor_(NULL)
{}

void
ScheduledRemoteSeedReader::onState() {
  if (notifier()->status() != Seed::INACTIVE && status() == READY) {
    statusIs(SCHEDULED);
    Activity * a = sendReactor_->notifier().ptr();
    IterationRank iteration = a->currentPhase() <= a->phase() ? a->currentIteration() : a->currentIteration().next();
    a->iterationIs(iteration);
    a->statusIs(Activity::scheduled);
  }
}

void
ScheduledRemoteSeedReader::phaseIs(PhaseRank p) {
  phase_ = p;
  sendReactor_->notifier()->phaseIs(p);
}

void
ScheduledRemoteSeedReader::SendReactor::onStatus() {
  switch (notifier()->status()) {
    case Activity::executing:
      parent()->statusIs(ScheduledRemoteSeedReader::BUSY);
      break;
    default:
      break;
  }
}

} // end namespace Pita
