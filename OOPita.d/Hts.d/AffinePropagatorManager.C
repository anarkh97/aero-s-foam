#include "AffinePropagatorManager.h"

#include "../AffineIntegratorPropagator.h"

#include <cassert>

namespace Pita { namespace Hts {

AffinePropagatorManager::AffinePropagatorManager(BasisCollector * collector,
                                                 GenFineIntegratorManager<AffineGenAlphaIntegrator> * integratorMgr,
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

AffineDynamPropagator *
AffinePropagatorManager::instance(const HalfSliceId & id) const {
  DynamPropagator * original = const_cast<DynamPropagator *>(collector_->source(id));
  if (original) {
    AffineDynamPropagator * downcasted = dynamic_cast<AffineDynamPropagator *>(original);
    assert(downcasted);
    return downcasted;
  }

  return NULL;
}

size_t
AffinePropagatorManager::instanceCount() const { 
  return collector_->sourceCount();
}

AffineDynamPropagator *
AffinePropagatorManager::instanceNew(const HalfSliceId & id) {
  // Create propagator 
  AffineGenAlphaIntegrator::Ptr integrator = integratorMgr_->fineIntegrator(id.direction());
  AffineIntegratorPropagator::Ptr newPropagator = AffineIntegratorPropagator::New(integrator.ptr());

  // Set up propagator
  Seconds halfCoarseTimeStep = fineTimeStep_ * halfSliceRatio_.value();
  HalfSliceRank initialSeedRank = (id.direction() == FORWARD) ? id.rank() : id.rank() + HalfSliceCount(1);
  Seconds sliceInitialTime = initialTime_ + halfCoarseTimeStep * initialSeedRank.value();

  newPropagator->initialTimeIs(sliceInitialTime);
  newPropagator->timeStepCountIs(halfSliceRatio_);
  
  // Attach BasisCollector Reactor 
  collector_->sourceIs(id, newPropagator.ptr());
 
  // Attach PostProcessing Reactor
  if (postProcessingMgr_) {
    this->postProcessingMgr_->outputFileSetIs(newPropagator.ptr(), PostProcessor::FileSetId(id.rank().value()));
  }

  return newPropagator.ptr();
}

void
AffinePropagatorManager::instanceDel(const HalfSliceId & id) {
  collector_->sourceIs(id, NULL);
}

} /* end namespace Hts */ } /* end namespace Pita */
