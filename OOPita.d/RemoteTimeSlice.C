#include "RemoteTimeSlice.h"

namespace Pita {

// RemoteTimeSlice implementation

  RemoteTimeSlice::RemoteTimeSlice(SliceRank rank, Status status, CpuRank owningCpu, Manager * manager) :
  TimeSliceImpl(rank, status),
  owningCpu_(owningCpu),
  manager_(manager),
  correctionReactor_( new CorrectionReactor( activityManagerInstance()->activityNew( String("RemoteCorrectionReactor ") + toString(rank.value()) ).ptr(), this ) )
{}

void
RemoteTimeSlice::onSeed() {
  if (status() == active) {
    manager()->outgoingCorrectionIs(initialInterface()->seed());
  }
}

void
RemoteTimeSlice::statusIs(Status s) {
  if (status() != s) {
    if (s == TimeSlice::stopped || s == TimeSlice::converged)
      finalInterface()->statusIs(TimeSlice::FinalInterface::final);
    setStatus(s);
  }
}

// RemoteTimeSlice::Manager implementation

RemoteTimeSlice::Manager::Manager(CommManager * commManager, size_t vectorSize) :
  vectorSize_(vectorSize),
  commManager_(commManager),
  commManagerReactor_(new CommManagerReactor(commManager, this)),
  commStatus_(available)
{}

RemoteTimeSlice::Ptr
RemoteTimeSlice::Manager::timeSlice(SliceRank rank) const {
  SliceMap::const_iterator it = slice_.find(rank);
  return (it != slice_.end()) ? it->second : NULL;
}

RemoteTimeSlice::Ptr
RemoteTimeSlice::Manager::timeSliceNew(SliceRank rank, CpuRank owningCpu) {
  SliceMap::iterator it = slice_.lower_bound(rank);
  if (it != slice_.end() && it->first == rank)
    throw Fwk::NameInUseException();
  RemoteTimeSlice::Ptr newSlice = RemoteTimeSlice::New(rank, owningCpu, this);
  slice_.insert(it, std::make_pair(rank, newSlice));
  return newSlice.ptr(); 
}

void
RemoteTimeSlice::Manager::timeSliceDel(SliceRank rank) {
  SliceMap::iterator it = slice_.find(rank);
  if (it != slice_.end())
    slice_.erase(it);
}

RemoteTimeSlice *
RemoteTimeSlice::Manager::getNextTimeSlice(const TimeSlice * ts) const {
  SliceMap::const_iterator it = ts ? slice_.upper_bound(ts->rank()) : slice_.begin(); 
  return (it != slice_.end()) ? it->second.ptr() : NULL;
}

void
RemoteTimeSlice::Manager::outgoingCorrectionIs(const DynamState & state) {
  commManager_->nextCorrectionIs(state);
  outgoingCorrection_ = state;
}

void
RemoteTimeSlice::Manager::commStatusIs(CommStatus s) {
  if (commStatus() != s) {
    commStatus_ = s;
    if (commStatus() == waiting) {
      commManagerReactor_->notifierIs(commManager_);
      commManager_->correctionVectorSizeIs(vectorSize());
      commStatus_ = available;
    }
  }
}

// RemoteTimeSlice::Manager::CommManagerReactor implementation

void
RemoteTimeSlice::Manager::CommManagerReactor::onReceivedCorrection() {
  parent_->incomingCorrection_ = notifier()->receivedCorrection();
}

// RemoteTimeSlice::CorrectionReactor implementation

RemoteTimeSlice::CorrectionReactor::CorrectionReactor(Activity * notifier, RemoteTimeSlice * slice) :
  TimeSlice::CorrectionReactor(notifier),
  slice_(slice)
{}

void
RemoteTimeSlice::CorrectionReactor::onStatus() {
  if (notifier()->status() == Activity::executing) {
    slice()->manager()->commStatusIs(RemoteTimeSlice::Manager::waiting);
    slice()->finalInterface()->nextSeedIs(slice()->manager()->incomingCorrection());
  }
}

} // end namespace Pita
