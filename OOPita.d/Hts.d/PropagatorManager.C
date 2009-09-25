#include "PropagatorManager.h"

namespace Pita { namespace Hts {

PropagatorManager::PropagatorManager(BasisCollector * collector,
                                                       FineIntegratorManager * integratorMgr,
                                                       PostProcessing::Manager * postProcessingMgr,
                                                       TimeStepCount halfSliceRatio,
                                                       Seconds initialTime) :
  collector_(collector),
  integratorMgr_(integratorMgr),
  postProcessingMgr_(postProcessingMgr),
  fineTimeStep_(integratorMgr->fineTimeStepSize()),
  halfSliceRatio_(halfSliceRatio),
  initialTime_(initialTime)
{}

IntegratorPropagator *
PropagatorManager::instance(const HalfSliceId & id) const {
  return collector_->source(id);
}

size_t
PropagatorManager::instanceCount() const { 
  return collector_->sourceCount();
}

IntegratorPropagator *
PropagatorManager::instanceNew(const HalfSliceId & id) {
  // Set up integrator
  DynamTimeIntegrator * integrator = integratorMgr_->fineIntegrator(id.direction());
 
  // Retrieve propagator 
  IntegratorPropagator * newPropagator = collector_->sourceNew(id, integrator);
 
  // Attach PostProcessor
  if (postProcessingMgr_) {
    this->postProcessingMgr_->outputFileSetIs(newPropagator, PostProcessor::FileSetId(id.rank().value()));
  }

  // Set up propagator
  Seconds halfCoarseTimeStep = fineTimeStep_ * halfSliceRatio_.value();
  HalfSliceRank initialSeedRank = (id.direction() == HalfTimeSlice::FORWARD) ? id.rank() : id.rank() + HalfSliceCount(1);
  Seconds sliceInitialTime = initialTime_ + halfCoarseTimeStep * initialSeedRank.value();

  newPropagator->initialTimeIs(sliceInitialTime);
  newPropagator->timeStepCountIs(halfSliceRatio_);

  return newPropagator;
}

void
PropagatorManager::instanceDel(const HalfSliceId & id) {
  collector_->sourceDel(id);
}

} /* end namespace Hts */ } /* end namespace Pita */
