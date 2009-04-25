#ifndef PITA_LOCALTIMESLICE_H
#define PITA_LOCALTIMESLICE_H

#include "Fwk.h"
#include "Types.h"
#include "TimeSliceImpl.h"
#include "DynamState.h"
#include "DynamStatePlainBasis.h"
#include "CommManager.h"
#include "DynamPropagator.h"
#include "Activity.h"

#include <map>

namespace Pita {

class LocalTimeSlice : public TimeSliceImpl {
public:
  typedef Fwk::Ptr<LocalTimeSlice> Ptr;
  typedef Fwk::Ptr<const LocalTimeSlice> PtrConst;
 
  Seconds initialTime() const { return initialTime_; }
  Seconds finalTime() const { return finalTime_; } 

  virtual void statusIs(Status s);
  
  virtual void onSeed();
  
  DynamStatePlainBasis::PtrConst initialProjectionBasis() const { return initialProjectionBasis_; }
  virtual void initialProjectionBasisInc(DynamStateBasis::PtrConst basis) = 0; 
  
  class Manager;
  
  class CorrectionReactor : public TimeSlice::CorrectionReactor {
  public:
    typedef Fwk::Ptr<CorrectionReactor> Ptr;
    typedef Fwk::Ptr<const CorrectionReactor> PtrConst;

    virtual void onStatus();
    
    virtual LocalTimeSlice * slice() const { return slice_; }
    void sliceIs(const LocalTimeSlice::Ptr & s) { slice_ = s.ptr(); }
    
    CorrectionReactor(Activity * notifier, LocalTimeSlice * slice);

  private:
    LocalTimeSlice * slice_;
  };

  friend class CorrectionReactor;
  
  virtual CorrectionReactor * correctionReactor() const { return correctionReactor_.ptr(); }
  
protected:
  LocalTimeSlice(SliceRank rank, Seconds initialTime, Seconds finalTime, Status status);

  virtual Manager * manager() const = 0;
  
  virtual DynamPropagator * seedUpdatePropagator() const = 0;
  DynamPropagator * localPropagator() const { return localPropagator_.ptr(); }
  DynamState previousSeed() const { return previousSeed_; }
  DynamState propagatedState() const { return propagatedState_; }

  void setLocalPropagator(DynamPropagator * newPropagator) { localPropagator_ = newPropagator; }
  void setPreviousSeed(const DynamState & newSeed) { previousSeed_ = newSeed; }
  void setPropagatedState(const DynamState & newState) { propagatedState_ = newState; }
  void setInitialProjectionBasis(const DynamStatePlainBasis * basis) { initialProjectionBasis_ = basis; }

  // Implementation function
  void propagateSeed();
  virtual void activateSlice() {} 
  virtual void optimizeLastIteration() {}
 
  // Implementation classes 
  class LocalPropagationReactor;
  friend class LocalPropagationReactor;
  class BasisUpdateReactor;
  friend class BasisUpdateReactor;
  
  LocalPropagationReactor * localPropagationReactor() const { return localPropagationReactor_.ptr(); }
  
  LocalTimeSlice(const LocalTimeSlice &); // No implementation
  LocalTimeSlice & operator=(const LocalTimeSlice &); // No implementation
  
private:
  Seconds initialTime_;
  Seconds finalTime_;
  DynamPropagator::Ptr localPropagator_;
  DynamState previousSeed_;
  DynamState propagatedState_;
  DynamStatePlainBasis::PtrConst initialProjectionBasis_;
  Fwk::Ptr<LocalPropagationReactor> localPropagationReactor_;
  Fwk::Ptr<CorrectionReactor> correctionReactor_;
  Fwk::Ptr<BasisUpdateReactor> basisUpdateReactor_;
};

class LocalTimeSlice::LocalPropagationReactor : public Activity::Notifiee {
public:
  typedef Fwk::Ptr<LocalPropagationReactor> Ptr;
  typedef Fwk::Ptr<const LocalPropagationReactor> PtrConst;

  virtual void onStatus();

  LocalPropagationReactor(Activity * notifier, LocalTimeSlice * slice);

private:
  LocalTimeSlice * slice_;
};

class LocalTimeSlice::BasisUpdateReactor : public Activity::Notifiee {
public:
  typedef Fwk::Ptr<BasisUpdateReactor> Ptr;
  typedef Fwk::Ptr<const BasisUpdateReactor> PtrConst;

  virtual void onStatus();

  BasisUpdateReactor(Activity * notifier, LocalTimeSlice * slice);

private:
  LocalTimeSlice * slice_;
};

class LocalTimeSlice::Manager : public TimeSlice::Manager {
public:
  typedef Fwk::Ptr<LocalTimeSlice::Manager> Ptr;
  typedef Fwk::Ptr<const LocalTimeSlice::Manager> PtrConst;

  LocalTimeSlice::Ptr timeSlice(SliceRank rank) const { return findTimeSlice(rank); }
  virtual SliceCount timeSliceCount() const = 0;
  
  LocalTimeSlice::Ptr timeSliceNew(SliceRank rank) { return createTimeSlice(rank); }
  virtual void timeSliceDel(SliceRank rank) = 0;
  
  CommManager::Ptr commManager() const { return commManager_; }
  void commManagerIs(CommManager::Ptr cm) { setCommManager(cm.ptr()); }
  
  DynamStateBasis::PtrConst incomingInitialStateBasis() const { return incomingInitialStateBasis_; }
  DynamStateBasis::PtrConst outgoingInitialStateBasis() const { return outgoingInitialStateBasis_; }
 
  void lastOutgoingInitialStateIs(const DynamState & state);
  void incomingInitialStateBasisIs(DynamStateBasis::Ptr basis); 

  class CommunicationReactor;  
  Fwk::Ptr<CommunicationReactor> communicationReactor() const { return communicationReactor_; }
 
  friend class TimeSliceIteratorConst<LocalTimeSlice>;
  TimeSliceIteratorConst<LocalTimeSlice> timeSliceIteratorConst() const {
    return TimeSliceIteratorConst<LocalTimeSlice>(this);
  }

  friend class LocalTimeSlice; 
  
protected:
  explicit Manager(CommManager * commManager);
  
  virtual LocalTimeSlice * createTimeSlice(SliceRank rank) = 0;
  virtual LocalTimeSlice * findTimeSlice(SliceRank rank) const = 0;

  virtual LocalTimeSlice * getNextTimeSlice(const TimeSlice *) const = 0;
  
  void setIncomingInitialStateBasis(DynamStateBasis * basis) { incomingInitialStateBasis_ = basis; }
  void setOutgoingInitialStateBasis(DynamStatePlainBasis * basis) { outgoingInitialStateBasis_ = basis; }
  void setCommManager(CommManager * cm) { commManager_ = cm; }
 
  class CommManagerReactor;
  friend class CommunicationReactor;

private:
  DynamStateBasis::Ptr incomingInitialStateBasis_;
  DynamStatePlainBasis::Ptr outgoingInitialStateBasis_;
  CommManager::Ptr commManager_;
  Fwk::Ptr<CommManagerReactor> commManagerReactor_;
  Fwk::Ptr<CommunicationReactor> communicationReactor_;
};

class LocalTimeSlice::Manager::CommManagerReactor : public CommManager::Notifiee {
public:
  typedef Fwk::Ptr<CommManagerReactor> Ptr;
  typedef Fwk::Ptr<const CommManagerReactor> PtrConst;

  virtual void onReceivedBroadcastBasis();

  CommManagerReactor(CommManager * notifier, LocalTimeSlice::Manager * parent);

private:
  LocalTimeSlice::Manager * parent_;
};


class LocalTimeSlice::Manager::CommunicationReactor : public Activity::Notifiee {
public:
  typedef Fwk::Ptr<CommunicationReactor> Ptr;
  typedef Fwk::Ptr<const CommunicationReactor> PtrConst;

  virtual void onStatus();
  
  CommunicationReactor(Activity * notifier, LocalTimeSlice::Manager * parent);
  
private:
  LocalTimeSlice::Manager * parent_;
};

} // end namespace Pita

#endif /* PITA_LOCALTIMESLICE_H */
