#ifndef PITA_HTS_PROPAGATOR_MANAGER
#define PITA_HTS_PROPAGATOR_MANAGER

#include "Fwk.h"
#include "Types.h"

#include "../IntegratorPropagator.h"
#include "BasisCollector.h"
#include "FineIntegratorManager.h"
#include "../PostProcessingManager.h"

namespace Pita { namespace Hts {

// DynamPropagator::Manager implementation for Pita HalfSlice problems
// Attach the necessary Reactors to the produced propagators
class PropagatorManager : public Fwk::GenManagerInterface<DynamPropagator*, HalfSliceId> {
public:
  EXPORT_PTRINTERFACE_TYPES(PropagatorManager);
 
  virtual IntegratorPropagator * instance(const HalfSliceId & id) const;
  virtual size_t instanceCount() const;

  virtual IntegratorPropagator * instanceNew(const HalfSliceId & id);
  virtual void instanceDel(const HalfSliceId & id);

//protected:
  PropagatorManager(BasisCollector * collector,
                    FineIntegratorManager * integratorMgr,
                    PostProcessing::Manager * postProcessingMgr,
                    TimeStepCount halfSliceRatio,
                    Seconds initialTime);

private:
  BasisCollector * collector_;
  FineIntegratorManager::Ptr integratorMgr_;
  PostProcessing::Manager::Ptr postProcessingMgr_;
  Seconds fineTimeStep_;
  TimeStepCount halfSliceRatio_;
  Seconds initialTime_;

  DISALLOW_COPY_AND_ASSIGN(PropagatorManager);
};

} /* end namespace Hts */ } /* end namespace Pita */

#endif /* PITA_HTS_PROPAGATOR_MANAGER */
