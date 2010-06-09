#ifndef PITA_COMMMANAGER_H
#define PITA_COMMMANAGER_H

#include "Fwk.h"
#include "Types.h"
#include "TimeSliceMapping.h"
#include "DynamStateBasis.h"

namespace Pita {

class CommManager : public Fwk::PtrInterface<CommManager> {
public:
  typedef Fwk::Ptr<CommManager> Ptr;
  typedef Fwk::Ptr<const CommManager> PtrConst;

  enum Status {
    available = 0,
    broadcasting,
    receiving,
    sending
  };

  Status status() const { return status_; }

  CpuCount availableCpus() const { return sliceMap_->availableCpus(); }
  CpuRank localCpuRank() const { return localCpuRank_; }
  CpuRank nextCpuRank() const { return nextCpuRank_; }
  CpuRank previousCpuRank() const { return previousCpuRank_; }
  
  SliceCount activeSlices() const { return sliceMap_->activeSlices(); }
 
  // Global Communication 
  DynamStateBasis::Ptr localBroadcastBasis() const { return localBroadcastBasis_; }
  DynamStateBasis::Ptr receivedBroadcastBasis() const { return receivedBroadcastBasis_; }
  virtual void localBroadcastBasisIs(DynamStateBasis::Ptr basis) = 0; 

  //Local Communication
  size_t correctionVectorSize() const { return receivedCorrection_.vectorSize(); }
  virtual void correctionVectorSizeIs(size_t s) = 0;
  DynamState receivedCorrection() const { return receivedCorrection_; }
  DynamState nextCorrection() const { return nextCorrection_; }
  virtual void nextCorrectionIs(const DynamState & state) = 0;
  SliceCount convergedSliceCount() const { return convergedSliceCount_; }
  
  class Notifiee : public Fwk::BaseNotifiee<CommManager> {
  public:
    typedef Fwk::Ptr<Notifiee> Ptr;
    typedef Fwk::Ptr<const Notifiee> PtrConst;

    virtual void onStatus() {}
    virtual void onReceivedBroadcastBasis() {}
    virtual void onReceivedCorrection() {}

  protected:
    explicit Notifiee(CommManager * notifier) : Fwk::BaseNotifiee<CommManager>(notifier) {} 
  };
  
  Notifiee::Ptr lastNotifiee() const { return notifiee_; }
  virtual void lastNotifieeIs(Notifiee * notifiee) { if (notifiee != notifiee_) setLastNotifiee(notifiee); }
 
protected:
  CommManager(TimeSliceMapping * sliceMap, CpuRank myCpuRank, CpuRank nextCpuRank, CpuRank previousCpuRank);

  const TimeSliceMapping * sliceMap() const { return sliceMap_.ptr(); }
  TimeSliceMapping * sliceMap() { return sliceMap_.ptr(); }

  void setStatus(Status s) { status_ = s; if (notifiee_) notifiee_->onStatus(); }
  void setLocalBroadcastBasis(DynamStateBasis * basis) { localBroadcastBasis_ = basis; } 
  void setReceivedBroadcastBasis(DynamStateBasis * basis) { receivedBroadcastBasis_ = basis; if (notifiee_) notifiee_->onReceivedBroadcastBasis(); }
  void setNextCorrection(const DynamState & state) { nextCorrection_ = state; }
  void setReceivedCorrection(const DynamState & state) { receivedCorrection_ = state; if (notifiee_) notifiee_->onReceivedCorrection(); }
  void setLastNotifiee(Notifiee * notifiee) { notifiee_ = notifiee; }
  void setConvergedSliceCount(SliceCount s) { convergedSliceCount_ = s; }
  
private:
  Status status_;
  TimeSliceMapping::Ptr sliceMap_;
  CpuRank localCpuRank_;
  CpuRank nextCpuRank_;
  CpuRank previousCpuRank_;
  Notifiee * notifiee_;
  DynamStateBasis::Ptr localBroadcastBasis_;
  DynamStateBasis::Ptr receivedBroadcastBasis_;
  DynamState receivedCorrection_;
  DynamState nextCorrection_;
  SliceCount convergedSliceCount_;
};
  
} // end namespace Pita

#endif /* PITA_COMMMANAGER_H */
