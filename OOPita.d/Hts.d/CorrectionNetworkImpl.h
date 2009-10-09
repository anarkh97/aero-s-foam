#ifndef PITA_HTS_CORRECTIONNETWORKIMPL_H
#define PITA_HTS_CORRECTIONNETWORKIMPL_H

#include "Fwk.h"
#include "Types.h"

#include "CorrectionNetwork.h"

#include "SliceMapping.h"

#include "../DynamOps.h"

#include "BasisCollectorImpl.h"

#include "../JumpProjectorImpl.h"
#include "../UpdatedSeedAssemblerImpl.h"
#include "ReducedFullTimeSliceImpl.h"

#include "../DynamStateReductorImpl.h"
#include "../DynamStateReconstructorImpl.h"

#include "../PivotedCholeskySolver.h"
#include "../NearSymmetricSolver.h"
#include <Math.d/FullSquareMatrix.h>

#include "../DynamStatePlainBasis.h"
#include "../SimpleBuffer.h"

#include <list>
#include <map>

class Communicator;

namespace Pita { namespace Hts {

class CorrectionNetworkImpl : public CorrectionNetwork {
public:
  EXPORT_PTRINTERFACE_TYPES(CorrectionNetworkImpl);

  virtual size_t reducedBasisSize() const { return metricBasis_->stateCount(); }

  virtual BasisCollector * collector() const { return collector_.ptr(); }
  
  virtual JumpProjectorImpl::Manager * jumpProjectorMgr() const { return jumpProjectorMgr_.ptr(); }
  virtual UpdatedSeedAssemblerImpl::Manager * updatedSeedAssemblerMgr() const { return updatedSeedAssemblerMgr_.ptr() ;} 
  virtual ReducedFullTimeSliceImpl::Manager * fullTimeSliceMgr() const { return fullTimeSliceMgr_.ptr(); }
  
  virtual DynamStateReductor::Manager * reductorMgr() const { return reductorMgr_.ptr(); }
  virtual DynamStateReconstructor::Manager * reconstructorMgr() const { return reconstructorMgr_.ptr(); }

  static Ptr New(size_t vSize, Communicator * timeComm, CpuRank myCpu,
                 const SliceMapping * mapping,
                 BasisCollectorImpl * collector,
                 const DynamOps * metric,
                 double projTol) {
    return new CorrectionNetworkImpl(vSize, timeComm, myCpu, mapping, collector, metric, projTol);
  }


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
    size_t stateCount() const;
    size_t stateCount(HalfTimeSlice::Direction d) const;
    size_t stateCount(CpuRank c) const;
    size_t stateCount(CpuRank c, HalfTimeSlice::Direction d) const;
    
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

  // HACK
  void buildProjection();

protected:
  CorrectionNetworkImpl(size_t vSize,
                        Communicator * timeComm,
                        CpuRank myCpu,
                        const SliceMapping * mapping,
                        BasisCollectorImpl * collector,
                        const DynamOps * metric,
                        double projectionTolerance);

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

  typedef std::map<int, DynamState> LocalBasis;
  LocalBasis localBasis_;
  DynamStatePlainBasis::Ptr metricBasis_;
  DynamStatePlainBasis::Ptr finalBasis_;

  FullSquareMatrix normalMatrix_;
  FullSquareMatrix reprojectionMatrix_;

  NearSymmetricSolver::Ptr solver_;
  
  BasisCollectorImpl::Ptr collector_;
  
  JumpProjectorImpl::Manager::Ptr jumpProjectorMgr_;
  UpdatedSeedAssemblerImpl::Manager::Ptr updatedSeedAssemblerMgr_;
  ReducedFullTimeSliceImpl::Manager::Ptr fullTimeSliceMgr_;

  DynamStateReductorImpl::Manager::Ptr reductorMgr_;
  DynamStateReconstructorImpl::Manager::Ptr reconstructorMgr_;

  typedef std::list<GlobalExchangeNumbering::Ptr> NumberingList;
  NumberingList globalExchangeNumbering_;
};

} /* end namespace Hts */ } /* end namespace Pita */

#endif /* PITA_HTS_CORRECTIONNETWORK_H */
