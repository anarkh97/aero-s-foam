#include "SeedUpdater.h"

#include <cassert>

namespace Pita {

void
SeedUpdater::iterationIs(IterationRank ir) {
  assert(propagatedSeed()->iteration() == ir);
  updateSeed();
  setIteration(ir);
}

void
SeedUpdater::updateSeed() {
  assert(correction()->iteration() == propagatedSeed()->iteration() || correction()->status() == Seed::INACTIVE);
  assert(correction()->status() != Seed::INACTIVE || propagatedSeed()->status() == Seed::CONVERGED); 

  DynamState newSeed = propagatedSeed()->state();
  if (propagatedSeed()->status() == Seed::ACTIVE) {
    newSeed += correction()->state();
  }

  updatedSeed()->stateIs(newSeed);
  updatedSeed()->statusIs(propagatedSeed()->status());
  updatedSeed()->iterationIs(propagatedSeed()->iteration());
}

size_t
SeedUpdater::vectorSize() const {
  return propagatedSeed() ? propagatedSeed()->state().vectorSize() : 0; 
}

} /* end namespace Pita */
