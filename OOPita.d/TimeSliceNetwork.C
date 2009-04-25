#include "TimeSliceNetwork.h"
#include "LocalTimeSlice.h"
#include "RemoteTimeSlice.h"
#include "DynamStatePlainBasis.h"

namespace Pita {

// Implementation helper classes
// Reactor to update the local slices

TimeSliceNetwork::ConvergenceReactor::ConvergenceReactor(Activity * notifier, TimeSliceNetwork * parent) :
  Activity::Notifiee(notifier),
  parent_(parent)
{}

void
TimeSliceNetwork::ConvergenceReactor::onStatus() {
  if (notifier()->status() == Activity::executing) {
    // TODO: allow for more than one slice converging simultaneously
    parent_->sliceMapping()->convergedSlicesInc();
    
    // Adjust Head Slice
    TimeSliceNetwork::SliceMap & sliceMap = parent_->slice_;
    SliceCount convergedSliceCount = parent_->sliceMapping()->convergedSlices();
    while (!sliceMap.empty() && (sliceMap.front()->rank().value() < convergedSliceCount.value())) {
      sliceMap.front()->statusIs(TimeSlice::converged);
      sliceMap.pop_front();
      // TODO Add inactive slices to sliceMap using sliceIt_
    }
  
    if (parent_->lastIteration() != notifier()->currentIteration()) {
      parent_->scheduleActivities();
      // Schedule next global communication
      Activity::Ptr activity = parent_->localSliceManager()->communicationReactor()->notifier();
      activity->iterationIs(activity->currentIteration());
      activity->statusIs(Activity::scheduled);
    } else {
      parent_->forceConvergence();
    }
  }
}


// TimeSliceNetwork implementation

TimeSliceNetwork::TimeSliceNetwork(CpuRank myCpu, TimeSliceMapping * mapping, LocalTimeSlice::Manager * localManager, RemoteTimeSlice::Manager * remoteManager, SeedInitializer * initializer) :
  myCpu_(myCpu),
  mapping_(mapping),
  sliceIt_(mapping_->slices(myCpu_)),
  localManager_(localManager),
  remoteManager_(remoteManager),
  initializer_(initializer),
  lastIteration_(0),
  convergenceReactor_( new ConvergenceReactor(activityManagerInstance()->activityNew("ConvergenceReactor").ptr(), this) )
{
  SliceRank firstActiveRank = mapping->firstActiveSlice();
  SliceRank firstInactiveRank = mapping->firstInactiveSlice();
  
  // Initial base
  DynamStatePlainBasis::Ptr initialProjectionBasis = DynamStatePlainBasis::New(initializer_->vectorSize());
  for (int i = 0; i < firstInactiveRank.value(); ++i) {
    initialProjectionBasis->lastStateIs(initializer_->initialSeed(SliceRank(i)));
  }

  // Initialize network
  LocalTimeSlice::Ptr newTimeSlice;
  
  for (; sliceIt_; ++sliceIt_) {
    SliceRank currentRank = *sliceIt_;
    if (currentRank >= firstInactiveRank)
      break;
    if (currentRank >= firstActiveRank) {
      newTimeSlice = localManager_->timeSliceNew(currentRank);
      if (currentRank > firstActiveRank) {
        SliceRank precedingRank = currentRank - SliceRank(1);
        TimeSlice::Ptr precedingSlice = localManager_->timeSlice(precedingRank);
        if (!precedingSlice) {
          precedingSlice = remoteManager_->timeSliceNew(precedingRank, mapping_->owningCpu(precedingRank));
          slice_.push_back(precedingSlice);
        }
        newTimeSlice->initialInterface()->notifierIs(precedingSlice->finalInterface());
        precedingSlice->finalInterface()->nextSeedIs(initializer_->initialSeed(currentRank));
        newTimeSlice->statusIs(TimeSlice::active);
      } else {
        newTimeSlice->initialInterface()->seedIs(initializer_->initialSeed(currentRank));
        newTimeSlice->statusIs(TimeSlice::lastIteration);
      }
      newTimeSlice->initialProjectionBasisInc(initialProjectionBasis);
      slice_.push_back(newTimeSlice);
    }
  }
  SliceRank nextRank = newTimeSlice->rank() + SliceRank(1);
  if (nextRank.value() < mapping->totalSlices().value()) {
    TimeSlice::Ptr followingSlice = remoteManager_->timeSliceNew(nextRank, mapping_->owningCpu(nextRank));
    followingSlice->initialInterface()->notifierIs(newTimeSlice->finalInterface());
    newTimeSlice->finalInterface()->nextSeedIs(initializer_->initialSeed(nextRank));
    followingSlice->statusIs(TimeSlice::active);
  }

  // Set up next correction for first active slice on cpu (whether local or remote)
  convergenceReactor_->notifier()->phaseIs(PhaseRank::convergence());
  convergenceReactor_->notifier()->iterationIs(convergenceReactor_->notifier()->currentIteration());
  scheduleActivities();
}

void
TimeSliceNetwork::scheduleActivities() {
  Activity::Ptr activity;
  
  // Schedule next correction
  if (!slice_.empty()) {
    activity = slice_.front()->correctionReactor()->notifier();
    activity->iterationIs(activity->currentIteration().next());
    activity->statusIs(Activity::scheduled);
  }
  
  // Schedule convergence check
  activity = convergenceReactor_->notifier();
  activity->iterationIs(activity->iteration().next());
  activity->statusIs(Activity::scheduled);
}

void
TimeSliceNetwork::forceConvergence() {
  // Force convergence after next fine grid computation
  TimeSliceNetwork::SliceMap::iterator end = slice_.end();
  for (TimeSliceNetwork::SliceMap::iterator it = slice_.begin(); it != end; ++it) {
    (*it)->statusIs(TimeSlice::lastIteration);
  }
}

} // end namespace Pita
