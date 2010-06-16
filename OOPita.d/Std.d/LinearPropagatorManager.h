#ifndef PITA_STD_LINEARPROPAGATORMANAGER_H
#define PITA_STD_LINEARPROPAGATORMANAGER_H

#include "Fwk.h"
#include "Types.h"

#include "../LinearDynamOps.h"
#include "../AffineIntegratorPropagator.h"
#include "../HomogeneousGenAlphaIntegrator.h"
#include "../PostProcessingManager.h"
#include "AffineBasisCollector.h"

namespace Pita { namespace Std {

class LinearPropagatorManager : public Fwk::GenManagerInterface<AffineDynamPropagator *, SliceRank> {
public:
  EXPORT_PTRINTERFACE_TYPES(LinearPropagatorManager);

  virtual AffineIntegratorPropagator * instance(const SliceRank & id) const;
  virtual size_t instanceCount() const;

  virtual AffineIntegratorPropagator * instanceNew(const SliceRank & id);
  virtual void instanceDel(const SliceRank & id);

  static Ptr New(AffineGenAlphaIntegrator * sharedIntegrator,
                 PostProcessing::Manager * postProcessingMgr,
                 AffineBasisCollector * collector,
                 TimeStepCount timeSliceRatio,
                 Seconds initialTime) {
    return new LinearPropagatorManager(sharedIntegrator, postProcessingMgr, collector, timeSliceRatio, initialTime);
  }

protected:
  LinearPropagatorManager(AffineGenAlphaIntegrator * sharedIntegrator,
                          PostProcessing::Manager * postProcessingMgr,
                          AffineBasisCollector * collector,
                          TimeStepCount timeSliceRatio,
                          Seconds initialTime);

  AffineIntegratorPropagator::Ptr createNewInstance(SliceRank rank);

private:
  AffineGenAlphaIntegrator::Ptr sharedIntegrator_;
  Seconds fineTimeStep_;
  TimeStepCount timeSliceRatio_;
  Seconds initialTime_;

  PostProcessing::Manager::Ptr postProcessingMgr_;
  AffineBasisCollector * collector_;
  
  typedef std::map<SliceRank, AffineIntegratorPropagator::Ptr> InstanceMap;
  InstanceMap instance_;
};

} /* end namespace Std */ } /* end namespace Pita */

#endif /* PITA_STD_LINEARPROPAGATORMANAGER_H */
