#ifndef PITA_REMOTETIMESLICE_H
#define PITA_REMOTETIMESLICE_H

#include "Fwk.h"
#include "Types.h"
#include "TimeSliceImpl.h"
#include "CommManager.h"
#include <map>

namespace Pita {

class RemoteTimeSlice : public TimeSliceImpl {
public:
  typedef Fwk::Ptr<RemoteTimeSlice> Ptr;
  typedef Fwk::Ptr<const RemoteTimeSlice> PtrConst;

  virtual void statusIs(Status s);

  CpuRank owningCpu() const { return owningCpu_; } 
  virtual void owningCpuIs(CpuRank cpu) { setOwningCpu(cpu); } 
 
  virtual void onSeed();
  
  class Manager;
 
  class CorrectionReactor : public TimeSlice::CorrectionReactor {
  public:
    typedef Fwk::Ptr<CorrectionReactor> Ptr;
    typedef Fwk::Ptr<const CorrectionReactor> PtrConst;

    virtual void onStatus();
    
    virtual RemoteTimeSlice * slice() const { return slice_; }
    void sliceIs(const RemoteTimeSlice::Ptr & s) { slice_ = s.ptr(); }
    
    CorrectionReactor(Activity * notifier, RemoteTimeSlice * slice);

  private:
    RemoteTimeSlice * slice_;
  };

  virtual CorrectionReactor * correctionReactor() const { return correctionReactor_.ptr(); }

  static RemoteTimeSlice::Ptr New(SliceRank rank, CpuRank owningCpu, RemoteTimeSlice::Manager * manager) { 
    return new RemoteTimeSlice(rank, TimeSlice::inactive, owningCpu, manager);
  }
  
protected:
  RemoteTimeSlice(SliceRank rank, Status status, CpuRank owningCpu, Manager * manager);
  
  Manager * manager() const { return manager_; }
  
  void setOwningCpu(CpuRank cpu) { owningCpu_ = cpu; } 
  
private:
  CpuRank owningCpu_;
  Manager * manager_;
  CorrectionReactor::Ptr correctionReactor_;
};


class RemoteTimeSlice::Manager : public Fwk::PtrInterface<RemoteTimeSlice::Manager> {
public:
  typedef Fwk::Ptr<RemoteTimeSlice::Manager> Ptr;
  typedef Fwk::Ptr<const RemoteTimeSlice::Manager> PtrConst;

  RemoteTimeSlice::Ptr timeSlice(SliceRank rank) const;
  SliceCount timeSliceCount() const { return SliceCount(slice_.size()); }
  
  RemoteTimeSlice::Ptr timeSliceNew(SliceRank rank, CpuRank owningCpu);
  void timeSliceDel(SliceRank rank);

  friend class TimeSliceIteratorConst<RemoteTimeSlice>;
  TimeSliceIteratorConst<RemoteTimeSlice> timeSliceIteratorConst() const {
    return TimeSliceIteratorConst<RemoteTimeSlice>(this);
  }

  DynamState outgoingCorrection() const { return outgoingCorrection_; }
  void outgoingCorrectionIs(const DynamState & state);

  enum CommStatus {
    available = 0,
    waiting
  };

  CommStatus commStatus() const { return commStatus_; }
  void commStatusIs(CommStatus s);

  DynamState incomingCorrection() const { return incomingCorrection_; }

  size_t vectorSize() const { return vectorSize_; }

  SliceCount convergedSliceCount() const { return convergedSliceCount_; }
  
  static RemoteTimeSlice::Manager::Ptr New(CommManager * commManager, size_t vectorSize) {
    return new Manager(commManager, vectorSize);
  }

protected:
  Manager(CommManager * commManager, size_t vectorSize);

  class CommManagerReactor;
  
  RemoteTimeSlice * getNextTimeSlice(const TimeSlice *) const;
 
private:
  typedef std::map<SliceRank, RemoteTimeSlice::Ptr> SliceMap;
  SliceMap slice_;
  size_t vectorSize_;
  CommManager::Ptr commManager_;
  Fwk::Ptr<CommManagerReactor> commManagerReactor_;
  CommStatus commStatus_;
  DynamState outgoingCorrection_;
  DynamState incomingCorrection_;
  SliceCount convergedSliceCount_;
};


class RemoteTimeSlice::Manager::CommManagerReactor : public CommManager::Notifiee {
public:
  typedef Fwk::Ptr<CommManagerReactor> Ptr;
  typedef Fwk::Ptr<const CommManagerReactor> PtrConst;

  virtual void onReceivedCorrection();

  CommManagerReactor(CommManager * notifier, RemoteTimeSlice::Manager * parent) :
    CommManager::Notifiee(notifier),
    parent_(parent)   
  {}

private:
  RemoteTimeSlice::Manager * parent_;
};

} // end namespace Pita

#endif /* PITA_REMOTETIMESLICE_H */
