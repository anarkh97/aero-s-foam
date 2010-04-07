#include "LocalTimeSlice.h"
#include "Activity.h"

#include <iostream>

namespace Pita {
  
// LocalTimeSlice implementation

LocalTimeSlice::LocalTimeSlice(SliceRank rank, Seconds initialTime, Seconds finalTime, Status status) :
  TimeSliceImpl(rank, status),
  initialTime_(initialTime),
  finalTime_(finalTime),
  localPropagator_(NULL),
  propagatedState_(),
  initialProjectionBasis_(NULL),
  localPropagationReactor_( new LocalPropagationReactor( activityManagerInstance()->activityNew( String("LocalPropagatorReactor ") + toString(rank.value()) ).ptr(), this ) ),
  correctionReactor_( new CorrectionReactor( activityManagerInstance()->activityNew( String("LocalCorrectionReactor ") + toString(rank.value()) ).ptr(), this ) ),
  basisUpdateReactor_( new BasisUpdateReactor( activityManagerInstance()->activityNew( String("BasisUpdateReactor ") + toString(rank.value()) ).ptr(), this ) )
{
  localPropagationReactor_->notifier()->phaseIs(PhaseRank::fineGrid());
  correctionReactor_->notifier()->phaseIs(PhaseRank::correction());
  basisUpdateReactor_->notifier()->phaseIs(PhaseRank::basisUpdate());
}

void
LocalTimeSlice::statusIs(Status s) {
  if (status() == s)
    return;

  bool currentlyActive = (status() == TimeSlice::active) || (status() == TimeSlice::lastIteration);
  bool becomingActive = (s == TimeSlice::active) || (s == TimeSlice::lastIteration);

  if ((!currentlyActive) && becomingActive) {
    activateSlice();
    basisUpdateReactor_->notifier()->iterationIs(basisUpdateReactor_->notifier()->currentIteration().next());
    basisUpdateReactor_->notifier()->iterationIs(localPropagationReactor_->notifier()->currentIteration());
    localPropagationReactor_->notifier()->statusIs(Activity::scheduled);
  }

  if (currentlyActive && (!becomingActive)) {
    finalInterface()->statusIs(TimeSlice::FinalInterface::final);
  }

  if (s == TimeSlice::lastIteration) {
    // No correction necessary
    optimizeLastIteration();
  }

  setStatus(s);
}

void
LocalTimeSlice::onSeed() {
  DynamState improvedSeed = initialInterface()->seed();
  if (status() == TimeSlice::active) {
    /* Correction */
    DynamState update = improvedSeed - previousSeed();
    setPreviousSeed(improvedSeed);
    seedUpdatePropagator()->initialStateIs(update);
    DynamState improvedNextSeed = propagatedState() + seedUpdatePropagator()->finalState();
    finalInterface()->nextSeedIs(improvedNextSeed);

    /* Schedule basis update and local fine propagation */
    basisUpdateReactor_->notifier()->statusIs(Activity::scheduled);
    localPropagationReactor()->notifier()->statusIs(Activity::scheduled);
    
    if (initialInterface()->status() == TimeSlice::FinalInterface::final) {
      statusIs(TimeSlice::lastIteration);
    }

  } else {
    setPreviousSeed(improvedSeed);
    if (status() == TimeSlice::lastIteration) {
      statusIs(TimeSlice::converged);
    }
  }
}

void
LocalTimeSlice::propagateSeed() {
  localPropagator()->initialStateIs(previousSeed());
  setPropagatedState(localPropagator()->finalState());
}

// Reactor implementations

LocalTimeSlice::LocalPropagationReactor::LocalPropagationReactor(Activity * notifier, LocalTimeSlice * slice) :
  Activity::Notifiee(notifier),
  slice_(slice)
{}
  
void
LocalTimeSlice::LocalPropagationReactor::onStatus() {
  if (slice_) {
    switch (notifier()->status()) {
      case Activity::executing:
        slice_->propagateSeed();
        break;
      case Activity::free:
         notifier()->iterationIs(notifier()->iteration().next()); // New iteration
         break;
      default:
         break;
    }
  }
}


LocalTimeSlice::CorrectionReactor::CorrectionReactor(Activity * notifier, LocalTimeSlice * slice) :
  TimeSlice::CorrectionReactor(notifier),
  slice_(slice)
{}

void
LocalTimeSlice::CorrectionReactor::onStatus() {
  if (slice_ && notifier()->status() == Activity::executing) {
    slice_->statusIs(TimeSlice::converged);
    slice_->finalInterface()->nextSeedIs(slice_->propagatedState());
  }
}

LocalTimeSlice::BasisUpdateReactor::BasisUpdateReactor(Activity * notifier, LocalTimeSlice * slice) :
  Activity::Notifiee(notifier),
  slice_(slice)
{}

void
LocalTimeSlice::BasisUpdateReactor::onStatus() {
  if (slice_) {
    switch (notifier()->status()) {
      case Activity::executing:
        slice_->manager()->lastOutgoingInitialStateIs(slice_->initialInterface()->seed());
        break;
      case Activity::free:
        notifier()->iterationIs(notifier()->iteration().next()); // New iteration
        break;
      default:
        break;
    }
  }
}


// Manager implementation

LocalTimeSlice::Manager::Manager(CommManager * commManager) :
  incomingInitialStateBasis_(NULL),
  outgoingInitialStateBasis_(NULL),
  commManager_(commManager),
  commManagerReactor_(new CommManagerReactor(commManager, this)),
  communicationReactor_(new CommunicationReactor(activityManagerInstance()->activityNew("CommunicationReactor").ptr(), this))
{
  communicationReactor_->notifier()->phaseIs(PhaseRank::globalCommunication());
  communicationReactor_->notifier()->iterationIs(communicationReactor_->notifier()->currentIteration());
}

void
LocalTimeSlice::Manager::lastOutgoingInitialStateIs(const DynamState & state) {
  outgoingInitialStateBasis_->lastStateIs(state);
}

void
LocalTimeSlice::Manager::incomingInitialStateBasisIs(DynamStateBasis::Ptr basis) {
  outgoingInitialStateBasis_->stateBasisDel();
  for (TimeSliceIteratorConst<LocalTimeSlice> it = timeSliceIteratorConst(); it; ++it) {
    it->initialProjectionBasisInc(basis);
  }
  setIncomingInitialStateBasis(basis.ptr());
} 

// Manager::CommunicationReactor implementation

LocalTimeSlice::Manager::CommunicationReactor::CommunicationReactor(Activity * n, LocalTimeSlice::Manager * parent) :
  Activity::Notifiee(n),
  parent_(parent)
{}

void
LocalTimeSlice::Manager::CommunicationReactor::onStatus() {
  if (notifier()->status() == Activity::executing) {
    parent_->commManagerReactor_->notifierIs(parent_->commManager_); 
    parent_->commManager_->localBroadcastBasisIs(parent_->outgoingInitialStateBasis_);
  }
}

// Manager::CommManagerReactor Implementation

LocalTimeSlice::Manager::CommManagerReactor::CommManagerReactor(CommManager * notifier, LocalTimeSlice::Manager * parent) :
  CommManager::Notifiee(notifier),
  parent_(parent)   
{}

void
LocalTimeSlice::Manager::CommManagerReactor::onReceivedBroadcastBasis() {
  //log() << "Received a global basis update, size = " << notifier()->receivedBroadcastBasis()->stateCount() << "\n";
  parent_->incomingInitialStateBasisIs(notifier()->receivedBroadcastBasis().ptr());
}

} // end namespace Pita
