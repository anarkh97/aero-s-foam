#include "RemoteCoarseCorrectionServer.h"

namespace Pita { namespace Hts {

RemoteCoarseCorrectionServer::RemoteCoarseCorrectionServer(
    RemoteDynamPropagatorServer * server,
    const SliceMapping * mapping,
    PhaseRank correctionPhase) :
  status_(IDLE),
  server_(server),
  mapping_(mapping),
  correctionPhase_(correctionPhase),
  correctionReactor_(NULL)
{
  Activity::Ptr activity = activityManagerInstance()->activityNew("RemoteCoarseCorrection");
  correctionReactor_ = new CorrectionReactor(activity.ptr(), this);
}

void
RemoteCoarseCorrectionServer::statusIs(Status s) {
  if (s != status()) {
    if (s == ACTIVE) {
      Activity * activity = correctionReactor_->notifier().ptr();
      activity->iterationIs(activity->currentIteration());
      activity->phaseIs(correctionPhase());
      activity->statusIs(Activity::scheduled);
    }
    status_ = s;
  }
}

// CorrectionReactor

RemoteCoarseCorrectionServer::CorrectionReactor::CorrectionReactor(
    Activity * notifier, RemoteCoarseCorrectionServer * parent) :
  Activity::Notifiee(notifier),
  parent_(parent)
{}

void
RemoteCoarseCorrectionServer::CorrectionReactor::onStatus() {
  if (notifier()->status() == Activity::executing) {
    const SliceMapping * mapping = parent_->mapping();
    int seedCount = mapping->activeSlices().value() / 2;

    for (int s = 0; s < seedCount; ++s) {
      CpuRank clientCpu = mapping->hostCpu(SliceId(BACKWARD_HALF_SLICE, HalfSliceRank(s * 2 + 1)));
      parent_->server()->initialStateNew(clientCpu);
    }
  }

  parent_->statusIs(ACTIVE);
}

} /* end namespace Hts */ } /* end namespace Pita */
