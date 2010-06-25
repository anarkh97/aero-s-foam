#include "NlSeedInitializerImpl.h"
#include <OOPita.d/PitaNonLinDynam.h>
#include <algorithm>

namespace Pita {

NlSeedInitializerImpl::NlSeedInitializerImpl(PitaNonLinDynamic * pbDesc,
                                         DynamTimeIntegrator * ti,
                                         TimeStepCount tsbs) : 
  SeedInitializer(pbDesc->solVecInfo(), SliceRank(pbDesc->getInitSeedCount())),
  probDesc_(pbDesc),
  lastFileSlice_(pbDesc->getInitSeedCount()),
  timeStepsBetweenSeeds_(tsbs),
  integrator_(ti),
  state_(DynamStatePlainBasis::New(vectorSize()))
{
  DynamState seed(vectorSize());
  pbDesc->getInitSeed(seed, lastFileSlice_.value() - 1);
  state_->lastStateIs(seed);
}

DynamState
NlSeedInitializerImpl::initialSeed(SliceRank rank) const {
  if (rank < lastFileSlice_) {
    DynamState seed(vectorSize());
    probDesc_->getInitSeed(seed, rank.value());
    return seed;
  }
  if (integrator_) {
    size_t targetIndex = static_cast<size_t>(rank.value())
                       - static_cast<size_t>(lastFileSlice_.value()) + 1;
    
    Seconds initialIntegrationTime = integrator_->timeStepSize() *
      ((lastFileSlice_.value() - 1) * (timeStepsBetweenSeeds().value()));

    integrator_->initialConditionIs(state_->state(state_->stateCount() - 1), initialIntegrationTime);

    for (size_t index = state_->stateCount() - 1; index < targetIndex; ++index) {
      integrator_->timeStepCountInc(this->timeStepsBetweenSeeds());
      state_->lastStateIs(integrator_->currentState()); 
    }
    return state_->state(targetIndex);
  }
  return DynamState(vectorSize());
}
  
} // end namespace Pita

