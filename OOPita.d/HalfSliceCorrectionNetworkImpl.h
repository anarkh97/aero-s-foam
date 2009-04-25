#ifndef PITA_HALFSLICECORRECTIONNETWORKIMPL_H
#define PITA_HALFSLICECORRECTIONNETWORKIMPL_H

#include "Fwk.h"
#include "Types.h"

#include "HalfSliceCorrectionNetwork.h"

#include "HalfSliceSchedule.h"
#include "SliceMapping.h"

#include "DynamOps.h"

#include "HalfSliceBasisCollectorImpl.h"

#include "DynamStateReductorImpl.h"
#include "DynamStateReconstructorImpl.h"

#include "Activity.h"

#include "DynamStatePlainBasis.h"
#include "SimpleBuffer.h"

class Communicator;

namespace Pita {

// TODO HACK
using namespace Hts;

class HalfSliceCorrectionNetworkImpl : public HalfSliceCorrectionNetwork {
public:
  EXPORT_PTRINTERFACE_TYPES(HalfSliceCorrectionNetworkImpl);

  enum Strategy {
    HOMOGENEOUS = 0,
    NON_HOMOGENEOUS
  };

  virtual HalfSliceBasisCollector * collector() const { return collector_.ptr(); }
  virtual DynamStateReductor::Manager * reductorMgr() const { return reductorMgr_.ptr(); }
  virtual DynamStateReconstructor::Manager * reconstructorMgr() const { return reconstructorMgr_.ptr(); }

  Strategy strategy() const { return strategy_; }

  PhaseRank correctionPhase() const { return correctionPhase_; }
  PhaseRank schedulingPhase() const { return schedulingPhase_; }

  static Ptr New(size_t vSize, Communicator * timeComm, CpuRank myCpu,
                 const HalfSliceSchedule * schedule,
                 const SliceMapping * mapping,
                 const DynamOps * metric,
                 Strategy strategy = HOMOGENEOUS) {
    return new HalfSliceCorrectionNetworkImpl(vSize, timeComm, myCpu, schedule, mapping, metric, strategy);
  }

  // Implementation classes
  class ProjectionBuildingReactor : public Activity::Notifiee {
  public:
    EXPORT_PTRINTERFACE_TYPES(ProjectionBuildingReactor);

    HalfSliceCorrectionNetworkImpl * parent() const { return parent_; }

    virtual void onStatus(); // Overriden

    ProjectionBuildingReactor(Activity * notifier, HalfSliceCorrectionNetworkImpl * parent) :
      Activity::Notifiee(notifier),
      parent_(parent)
    {}

  private:
    HalfSliceCorrectionNetworkImpl * parent_;
  };

  class SchedulingReactor : public Activity::Notifiee {
  public:
    EXPORT_PTRINTERFACE_TYPES(SchedulingReactor);

    HalfSliceCorrectionNetworkImpl * parent() const { return parent_; }

    virtual void onStatus(); // Overriden

    SchedulingReactor(Activity * notifier, HalfSliceCorrectionNetworkImpl * parent) :
      Activity::Notifiee(notifier),
      parent_(parent)
    {}

  private:
    HalfSliceCorrectionNetworkImpl * parent_;
  };

  class NonHomogeneousSchedulingReactor : public SchedulingReactor {
  public:
    EXPORT_PTRINTERFACE_TYPES(NonHomogeneousSchedulingReactor);

    virtual void onStatus(); // Overriden

    NonHomogeneousSchedulingReactor(Activity * notifier, HalfSliceCorrectionNetworkImpl * parent) :
      SchedulingReactor(notifier, parent)
    {}
  };

  class GlobalExchangeNumbering : public Fwk::PtrInterface<GlobalExchangeNumbering> {
  public:
    EXPORT_PTRINTERFACE_TYPES(GlobalExchangeNumbering);

    class IteratorConst {
    public:
      std::pair<HalfTimeSlice::Direction, int> operator*() const { return std::make_pair(it_->first.direction(), it_->second); }
      IteratorConst & operator++() { ++it_; return *this; }
      IteratorConst operator++(int) { IteratorConst tmp(*this); this->operator++(); return tmp; }
      operator bool() const { return it_ != endIt_; }
     
      IteratorConst() {
        it_ = endIt_;
      }

    protected:
      typedef std::map<HalfSliceId, size_t>::const_iterator ItImpl;

      explicit IteratorConst(ItImpl beginIt, ItImpl endIt) : 
        it_(beginIt),
        endIt_(endIt)
      {}

      friend class GlobalExchangeNumbering;

    private:
      ItImpl it_;
      ItImpl endIt_;
    };
   
    // State count
    size_t stateCount(HalfTimeSlice::Direction d = HalfTimeSlice::NO_DIRECTION) const;
    size_t stateCount(CpuRank c, HalfTimeSlice::Direction d = HalfTimeSlice::NO_DIRECTION) const;
    
    // Numbering for all states (both initial AND final)
    int globalIndex(const HalfSliceId & id) const;
    IteratorConst globalIndex() const;
    HalfSliceId stateId(int globalFullIndex) const;
    
    // Numbering for initial OR final states
    int globalHalfIndex(const HalfSliceId & id) const;
    IteratorConst globalHalfIndex(HalfTimeSlice::Direction d) const;
    HalfSliceId stateId(int globalHalfIndex, HalfTimeSlice::Direction d) const;

    explicit GlobalExchangeNumbering(const SliceMapping * m);

  private:
    typedef std::map<HalfSliceId, size_t> IndexMap;

    void initialize(const SliceMapping * m);

    std::vector<size_t> stateCount_;
    std::vector<size_t> initialStateCount_;
    std::vector<size_t> finalStateCount_;
    
    IndexMap globalIndex_;
    IndexMap initialGlobalIndex_;
    IndexMap finalGlobalIndex_;
    
    std::vector<HalfSliceId> stateId_;
    std::vector<HalfSliceId> initialStateId_;
    std::vector<HalfSliceId> finalStateId_;

    DISALLOW_COPY_AND_ASSIGN(GlobalExchangeNumbering);
  };

protected:
  HalfSliceCorrectionNetworkImpl(size_t vSize,
                                 Communicator * timeComm,
                                 CpuRank myCpu,
                                 const HalfSliceSchedule * schedule,
                                 const SliceMapping * mapping,
                                 const DynamOps * metric,
                                 Strategy strategy);


  void buildProjection();

  friend class ProjectionBuildingReactor;
  friend class SchedulingReactor;

private:
  size_t vectorSize_;
  
  CpuRank localCpu_;
  Communicator * timeCommunicator_;

  PhaseRank correctionPhase_;
  PhaseRank schedulingPhase_;

  SliceMapping::PtrConst mapping_;
 
  DynamOps::PtrConst metric_;
  
  SimpleBuffer<double> gBuffer_;
  SimpleBuffer<double> mBuffer_;
  SimpleBuffer<int> mpiParameters_;

  DynamStatePlainBasis::Ptr metricBasis_;
  DynamStatePlainBasis::Ptr finalBasis_;

  SymFullMatrix normalMatrix_;
  PivotedCholeskySolver::Ptr solver_;
  
  HalfSliceBasisCollectorImpl::Ptr collector_;
  DynamStateReductorImpl::Manager::Ptr reductorMgr_;
  DynamStateReconstructorImpl::Manager::Ptr reconstructorMgr_;

  ProjectionBuildingReactor::Ptr projectionBuildingReactor_;
  Strategy strategy_;
  SchedulingReactor::Ptr schedulingReactor_;

  GlobalExchangeNumbering::Ptr globalExchangeNumbering_;
};

} // end namespace Pita

#endif /* PITA_HALFSLICECORRECTIONNETWORK_H */
