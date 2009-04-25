#include "IntegratorSeedInitializer.h"

namespace Pita {

IntegratorSeedInitializer::IntegratorSeedInitializer(DynamTimeIntegrator * i,
                                                     TimeStepCount tsbs) :
  SeedInitializer(i->vectorSize(), SliceRank(0)),
  integrator_(i),
  timeStepsBetweenSeeds_(tsbs),
  seed_(DynamStatePlainBasis::New(i->vectorSize()))
{
  seed_->lastStateIs(i->currentState());
}

DynamState
IntegratorSeedInitializer::initialSeed(SliceRank rank) const {
  if (rank > this->lastSlice()) {
    
    int targetIndex = rank.value();
    IntegratorSeedInitializer * mutable_this = const_cast<IntegratorSeedInitializer *>(this);

    for (int index = this->lastSlice().value(); index < targetIndex; ++index) {
      mutable_this->integrator_->timeStepCountInc(this->timeStepsBetweenSeeds());
      mutable_this->seed_->lastStateIs(integrator_->currentState()); 
    }

    mutable_this->setSlices(rank);
  }

  return seed_->state(rank.value());
}

} // end namespace Pita
