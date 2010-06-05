#ifndef PITA_HTS_LINEARPROJECTIONNETWORKIMPL_H
#define PITA_HTS_LINEARPROJECTIONNETWORKIMPL_H

#include "Fwk.h"
#include "Types.h"

#include "ProjectionNetwork.h"

#include "SliceMapping.h"
#include "../DynamOps.h"
#include "BasisCollectorImpl.h"
#include "../DynamStatePlainBasis.h"

#include "../SimpleBuffer.h"
#include <list>
#include <map>

class Communicator;

namespace Pita { namespace Hts {

class LinearProjectionNetworkImpl : public ProjectionNetwork {
public:
  EXPORT_PTRINTERFACE_TYPES(LinearProjectionNetworkImpl);

  virtual size_t reducedBasisSize() const { return metricBasis_->stateCount(); }

  virtual BasisCollector * collector() const { return collector_.ptr(); }
  
  virtual const DynamStateBasis * projectionBasis() const { return metricBasis_.ptr(); }
  virtual const DynamStateBasis * propagatedBasis() const { return finalBasis_.ptr(); }
  virtual const RankDeficientSolver * normalMatrixSolver() const { return solver_.ptr(); }
  virtual const FullSquareMatrix * reprojectionMatrix() const { return &reprojectionMatrix_; }
 
  const DynamOps * metric() const { return metric_.ptr(); }

  static Ptr New(size_t vSize, Communicator * timeComm, CpuRank myCpu,
                 const SliceMapping * mapping,
                 BasisCollectorImpl * collector,
                 const DynamOps * metric,
                 RankDeficientSolver * solver) {
    return new LinearProjectionNetworkImpl(vSize, timeComm, myCpu, mapping, collector, metric, solver);
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
  void prepareProjection();
  void buildProjection();

protected:
  LinearProjectionNetworkImpl(size_t vSize,
                        Communicator * timeComm,
                        CpuRank myCpu,
                        const SliceMapping * mapping,
                        BasisCollectorImpl * collector,
                        const DynamOps * metric,
                        RankDeficientSolver * solver);

private:
  size_t vectorSize_;
  
  CpuRank localCpu_;
  Communicator * timeCommunicator_;

  SliceMapping::PtrConst mapping_;
 
  DynamOps::PtrConst metric_;
  
  SimpleBuffer<double> gBuffer_;
  SimpleBuffer<double> mBuffer_;
  SimpleBuffer<int> mpiParameters_;

  typedef std::map<int, DynamState> LocalBasis;
  LocalBasis localBasis_;
  
  DynamStatePlainBasis::Ptr metricBasis_;
  DynamStatePlainBasis::Ptr finalBasis_;
  
  DynamStatePlainBasis::Ptr originalMetricBasis_;
  DynamStatePlainBasis::Ptr originalFinalBasis_;
  
  FullSquareMatrix normalMatrix_;
  FullSquareMatrix transmissionMatrix_;
  FullSquareMatrix reprojectionMatrix_;
  RankDeficientSolver::Ptr solver_;
  
  BasisCollectorImpl::Ptr collector_;
  
  typedef std::list<GlobalExchangeNumbering::Ptr> NumberingList;
  NumberingList globalExchangeNumbering_;
};

} /* end namespace Hts */ } /* end namespace Pita */

#endif /* PITA_HTS_LINEARPROJECTIONNETWORK_H */
