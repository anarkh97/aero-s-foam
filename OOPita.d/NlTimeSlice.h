#ifndef PITA_NLTIMESLICE_H
#define PITA_NLTIMESLICE_H

#include "Fwk.h"
#include "Types.h"
#include "DynamState.h"
#include "LocalTimeSlice.h"
#include "NlDynamTimeIntegrator.h"
#include "ProjectorPropagator.h"
#include "NlDynamOps.h"

namespace Pita {

class PitaNonLinDynamic;
  
class NlTimeSlice : public LocalTimeSlice {
public:
  typedef Fwk::Ptr<NlTimeSlice> Ptr;
  typedef Fwk::Ptr<const NlTimeSlice> PtrConst;
 
  virtual void initialProjectionBasisInc(DynamStateBasis::PtrConst basis);
  
  class Manager : public LocalTimeSlice::Manager {
  public:
    typedef Fwk::Ptr<NlTimeSlice::Manager> Ptr;
    typedef Fwk::Ptr<const NlTimeSlice::Manager> PtrConst;

    NlTimeSlice::Ptr timeSlice(SliceRank rank) const { return findTimeSlice(rank); }
    virtual SliceCount timeSliceCount() const { return SliceCount(slice_.size()); }
    
    NlTimeSlice::Ptr timeSliceNew(SliceRank rank) { return createTimeSlice(rank); }
    virtual void timeSliceDel(SliceRank rank);
   
    friend class TimeSliceIteratorConst<NlTimeSlice>;
    TimeSliceIteratorConst<NlTimeSlice> timeSliceIteratorConst() {
      return TimeSliceIteratorConst<NlTimeSlice>(this);
    }
    
    static NlTimeSlice::Manager::Ptr New(NlDynamTimeIntegrator * integrator, CommManager * commManager) { 
      return new Manager(integrator, commManager);
    }
    
  protected:
    Manager(NlDynamTimeIntegrator * integrator, CommManager * commManager);
    virtual NlTimeSlice * createTimeSlice(SliceRank rank);
    virtual NlTimeSlice * findTimeSlice(SliceRank rank) const;
   
    virtual NlTimeSlice * getNextTimeSlice(const TimeSlice * ts) const;
    
    Seconds coarseTimeStep() const;
    TimeStepCount timeGridRatio() const;
    size_t vectorSize() const;
    
    PitaNonLinDynamic * probDesc() const { return probDesc_; }
    
  private:
    typedef std::map<SliceRank, NlTimeSlice::Ptr > SliceMap;
    SliceMap slice_; 

    PitaNonLinDynamic * probDesc_;
    NlDynamOps::Ptr dynamOps_; 
    NlDynamTimeIntegrator::Ptr timeIntegrator_;
  };
  
  static NlTimeSlice::Ptr New(SliceRank rank, Seconds t0, Seconds tf, Manager * manager, Status status = inactive) {
    return new NlTimeSlice(rank, t0, tf, manager, status);
  }
 
protected:
  NlTimeSlice(SliceRank rank, Seconds initialTime, Seconds finalTime, Manager * manager, Status status);

  Manager * manager() const { return manager_; }

  virtual ProjectorPropagator * seedUpdatePropagator() const { return seedUpdatePropagator_.ptr(); }
  void setSeedUpdatePropagator(ProjectorPropagator * propagator) { seedUpdatePropagator_ = propagator; }
 
  virtual void activateSlice();
  virtual void optimizeLastIteration();
  
private:
  Manager * manager_;

  ProjectorPropagator::Ptr seedUpdatePropagator_;
};

} // end namespace Pita

#endif /* PITA_NLTIMESLICE_H */
