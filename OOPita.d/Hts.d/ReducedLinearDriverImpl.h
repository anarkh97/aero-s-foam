#ifndef PITA_HTS_REDUCEDLINEARDRIVERIMPL_H
#define PITA_HTS_REDUCEDLINEARDRIVERIMPL_H

#include "Fwk.h"
#include "Types.h"

#include "../LinearDriver.h"

#include "../DynamState.h"

#include "SliceMapping.h"

#include "../PostProcessingManager.h"
#include "BasisCollectorImpl.h"
#include "LocalCorrectionTimeSlice.h"
#include "../SeedInitializer.h"

template <typename Scalar> class SingleDomainDynamic;
class Domain;
class SolverInfo;
class Communicator;
class SDDynamPostProcessor;

namespace Pita { namespace Hts {

class ReducedLinearDriverImpl : public LinearDriver {
public:
  EXPORT_PTRINTERFACE_TYPES(ReducedLinearDriverImpl);

  virtual void solve();

  SingleDomainDynamic<double> * probDesc() const { return probDesc_; }
  Domain * domain() const { return domain_; }
  SolverInfo * solverInfo() const { return solverInfo_; }
  Communicator * baseComm() const { return baseComm_; }

  // Independent from global state
  static ReducedLinearDriverImpl::Ptr New(SingleDomainDynamic<double> * pbDesc,
                                          Domain * domain,
                                          SolverInfo * solverInfo,
                                          Communicator * baseComm) {
    return new ReducedLinearDriverImpl(pbDesc, domain, solverInfo, baseComm);
  }

protected:
  ReducedLinearDriverImpl(SingleDomainDynamic<double> *, Domain *, SolverInfo *, Communicator *);

  void preprocess();
  void solveParallel(Communicator * timeComm, Communicator * coarseComm);
  void solveCoarse(Communicator * timeComm);
 
  DynamState initialSeed() const;
  
  HalfTimeSlice::Manager::Ptr buildHalfTimeSliceManager(
      GeneralizedAlphaParameter fineIntegrationParam,
      PostProcessing::Manager * postProcessingMgr,
      BasisCollector * collector) const;
  
  PostProcessing::Manager::Ptr buildPostProcessor(CpuRank localCpu) const;
  BasisCollectorImpl::Ptr buildBasisCollector() const;
  CorrectionTimeSlice::Manager::Ptr buildCoarseCorrection(Communicator * coarseComm) const;
  LinearGenAlphaIntegrator::Ptr buildCoarseIntegrator() const; 
  DynamPropagator::Ptr buildCoarsePropagator(bool local, Communicator * coarseComm) const;
  SeedInitializer::Ptr buildSeedInitializer(bool local, Communicator * timeComm) const;

private:
  /* Primary sources */
  SingleDomainDynamic<double> * probDesc_;
  Domain * domain_;
  SolverInfo * solverInfo_;
  Communicator * baseComm_;

  /* Space-domain */
  size_t vectorSize_;
  LinearDynamOps::Manager::Ptr dynamOpsMgr_;

  /* Time-domain */ 
  Seconds fineTimeStep_;
  TimeStepCount halfSliceRatio_;
  TimeStepCount sliceRatio_;
  Seconds coarseTimeStep_;
  Seconds initialTime_;
  Seconds finalTime_;

  /* Main options */
  bool noForce_;
  bool remoteCoarse_;
  
  /* Load balancing */ 
  SliceMapping::Ptr mapping_;

  /* Other parameters */
  IterationRank lastIteration_;
  double projectorTolerance_;
  double coarseRhoInfinity_;
};

} /* namespace Hts */ } /* end namespace Pita */

Pita::LinearDriver::Ptr linearPitaDriverNew(SingleDomainDynamic<double> * pbDesc);

#endif /* PITA_HTS_REDUCEDLINEARDRIVERIMPL_H */
