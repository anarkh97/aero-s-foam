#ifndef PITA_HTS_PROPAGATOR_MANAGER
#define PITA_HTS_PROPAGATOR_MANAGER

#include "Fwk.h"
#include "Types.h"

#include "../AffineDynamPropagator.h"

#include "FineIntegratorManager.h"
#include "../PostProcessingManager.h"
#include "BasisCollector.h"

#include "../HomogeneousGenAlphaIntegrator.h"

namespace Pita { namespace Hts {

// DynamPropagator::Manager implementation for Pita HalfSlice problems
// Attach the necessary Reactors to the produced propagators
class PropagatorManager : public Fwk::GenManagerInterface<AffineDynamPropagator*, HalfSliceId> {
public:
  EXPORT_PTRINTERFACE_TYPES(PropagatorManager);
 
  virtual AffineDynamPropagator * instance(const HalfSliceId & id) const;
  virtual size_t instanceCount() const;

  virtual AffineDynamPropagator * instanceNew(const HalfSliceId & id);
  virtual void instanceDel(const HalfSliceId & id);

  static Ptr New(BasisCollector * collector,
                 GenFineIntegratorManager<AffineGenAlphaIntegrator> * integratorMgr,
                 PostProcessing::Manager * postProcessingMgr,
                 TimeStepCount halfSliceRatio,
                 Seconds initialTime) {
    return new PropagatorManager(collector, integratorMgr, postProcessingMgr, halfSliceRatio, initialTime);
  }

protected:
  PropagatorManager(BasisCollector * collector,
                    GenFineIntegratorManager<AffineGenAlphaIntegrator> * integratorMgr,
                    PostProcessing::Manager * postProcessingMgr,
                    TimeStepCount halfSliceRatio,
                    Seconds initialTime);

private:
  BasisCollector * collector_;
  GenFineIntegratorManager<AffineGenAlphaIntegrator>::Ptr integratorMgr_;
  PostProcessing::Manager::Ptr postProcessingMgr_;
  Seconds fineTimeStep_;
  TimeStepCount halfSliceRatio_;
  Seconds initialTime_;

  DISALLOW_COPY_AND_ASSIGN(PropagatorManager);
};

} /* end namespace Hts */ } /* end namespace Pita */

#endif /* PITA_HTS_PROPAGATOR_MANAGER */
