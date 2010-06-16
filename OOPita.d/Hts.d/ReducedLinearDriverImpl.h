#ifndef PITA_HTS_REDUCEDLINEARDRIVERIMPL_H
#define PITA_HTS_REDUCEDLINEARDRIVERIMPL_H

#include "Fwk.h"
#include "Types.h"

#include "../LinearDriverImpl.h"

#include "../DynamState.h"

#include "SliceMapping.h"

#include "../LinearGenAlphaIntegrator.h"
#include "../PostProcessingManager.h"
#include "BasisCollectorImpl.h"
#include "../CorrectionPropagator.h"
#include "../SeedInitializer.h"

template <typename Scalar> class SingleDomainDynamic;
class GeoSource;
class Domain;
class SolverInfo;
class Communicator;
class SDDynamPostProcessor;

namespace Pita { namespace Hts {

class ReducedLinearDriverImpl : public LinearDriverImpl {
public:
  EXPORT_PTRINTERFACE_TYPES(ReducedLinearDriverImpl);

  virtual void solve();
  
  static ReducedLinearDriverImpl::Ptr New(SingleDomainDynamic<double> * pbDesc,
                                          GeoSource * geoSource,
                                          Domain * domain,
                                          SolverInfo * solverInfo,
                                          Communicator * baseComm) {
    return new ReducedLinearDriverImpl(pbDesc, geoSource, domain, solverInfo, baseComm);
  }

protected:
  ReducedLinearDriverImpl(SingleDomainDynamic<double> *, GeoSource *, Domain *, SolverInfo *, Communicator *);

  void preprocess();
  void solveParallel(Communicator * timeComm, Communicator * coarseComm);
  void solveCoarse(Communicator * timeComm);
 
  HalfTimeSlice::Manager::Ptr buildHalfTimeSliceManager(
      GeneralizedAlphaParameter fineIntegrationParam,
      PostProcessing::Manager * postProcessingMgr,
      BasisCollector * collector) const;
  
  PostProcessing::Manager::Ptr buildPostProcessor(CpuRank localCpu) const;
  BasisCollectorImpl::Ptr buildBasisCollector() const;
  CorrectionPropagator<DynamState>::Manager::Ptr buildCoarseCorrection(Communicator * coarseComm) const;
  LinearGenAlphaIntegrator::Ptr buildCoarseIntegrator() const; 
  DynamPropagator::Ptr buildCoarsePropagator(bool local, Communicator * coarseComm) const;
  SeedInitializer::Ptr buildSeedInitializer(bool local, Communicator * timeComm) const;

private:
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

Pita::LinearDriver::Ptr linearReversiblePitaDriverNew(SingleDomainDynamic<double> * pbDesc);

#endif /* PITA_HTS_REDUCEDLINEARDRIVERIMPL_H */
