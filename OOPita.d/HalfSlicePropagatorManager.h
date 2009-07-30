#ifndef PITA_HALFSLICEPROPAGATOR_MANAGER
#define PITA_HALFSLICEPROPAGATOR_MANAGER

#include "Fwk.h"
#include "Types.h"

#include "IntegratorPropagator.h"
#include "HalfSliceBasisCollector.h"
#include "FineIntegratorManager.h"
#include "PostProcessingManager.h"

namespace Pita {

// DynamPropagator::Manager implementation for Pita HalfSlice problems
// Attach the necessary Reactors to the produced propagators
class HalfSlicePropagatorManager : public Fwk::GenManagerInterface<DynamPropagator*, HalfSliceId> {
public:
  EXPORT_PTRINTERFACE_TYPES(HalfSlicePropagatorManager);
 
  virtual IntegratorPropagator * instance(const HalfSliceId & id) const;
  virtual size_t instanceCount() const;

  virtual IntegratorPropagator * instanceNew(const HalfSliceId & id);
  virtual void instanceDel(const HalfSliceId & id);

//protected:
  HalfSlicePropagatorManager(HalfSliceBasisCollector * collector,
                             FineIntegratorManager * integratorMgr,
                             PostProcessing::Manager * postProcessingMgr,
                             TimeStepCount halfSliceRatio,
                             Seconds initialTime);

private:
  HalfSliceBasisCollector * collector_;
  FineIntegratorManager::Ptr integratorMgr_;
  PostProcessing::Manager::Ptr postProcessingMgr_;
  Seconds fineTimeStep_;
  TimeStepCount halfSliceRatio_;
  Seconds initialTime_;

  DISALLOW_COPY_AND_ASSIGN(HalfSlicePropagatorManager);
};

} // end namespace Pita

#endif /* PITA_HALFSLICEPROPAGATOR_MANAGER */
