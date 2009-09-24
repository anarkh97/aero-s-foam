#include "RemoteHalfSliceSurrogate.h"

namespace Pita {

RemoteHalfSliceSurrogate::RemoteHalfSliceSurrogate(SliceMapping * mapping) :
  mapping_(mapping),
  iteration_(activityManagerInstance()->currentIteration())
{}

void
RemoteHalfSliceSurrogate::iterationIs(IterationRank i) {
  iteration_ = i;
  RecvRectorContainer & reactor = (i.value() % 2 == 0) ? evenReactor_ : oddReactor_;
  for (RecvRectorContainer::iterator it = reactor.begin(); it != reactor.end(); ++it) {
    if (it->first.rank() > mapping_->firstActiveSlice() && it->first.rank() <= mapping_->firstInactiveSlice()) {
      Activity * a = it->second->notifier().ptr();
      a->iterationIs(a->currentIteration());
      a->statusIs(Activity::scheduled);
    }
  }
}

inline
RemoteHalfSliceSurrogate::RecvRectorContainer &
RemoteHalfSliceSurrogate::getReactorContainer(const Hts::SeedId & key) {
  // HACK ! TODO
  if (key.type() == Hts::UNDEFINED_SEED) {
    return (key.rank().value() % 2 == 0) ? oddReactor_ : evenReactor_;
  }
  return (key.rank().value() % 2 == 0) ? evenReactor_ : oddReactor_;
}

inline
const RemoteHalfSliceSurrogate::RecvRectorContainer &
RemoteHalfSliceSurrogate::getReactorContainer(const Hts::SeedId & key) const {
  return const_cast<RemoteHalfSliceSurrogate *>(this)->getReactorContainer(key); 
}

void
RemoteHalfSliceSurrogate::instanceNew(const Hts::SeedId & key, ScheduledRemoteSeedWriter * writer) {
  Activity * activity = activityManagerInstance()->activityNew("RemoteWriter_" + toString(key.type()) + toString(key.rank())).ptr();
  activity->iterationIs(activityManagerInstance()->currentIteration());
  ReceiveReactor::Ptr reactor = new ReceiveReactor(activity, writer);
  RecvRectorContainer & reactorContainer = getReactorContainer(key);
  reactorContainer[key] = reactor;
  setReceiveReactor(writer, reactor.ptr());
}

void
RemoteHalfSliceSurrogate::instanceDel(const Hts::SeedId & key) {
  RecvRectorContainer & reactorContainer = getReactorContainer(key);
  reactorContainer.erase(key); 
}

} // end namespace Pita
