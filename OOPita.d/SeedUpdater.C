#include "SeedUpdater.h"

#include <cassert>

namespace Pita {

void
SeedUpdater::iterationIs(IterationRank ir) {
  assert(propagatedSeed()->iteration().next() == ir);
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
  updatedSeed()->iterationIs(propagatedSeed()->iteration().next());
}

size_t
SeedUpdater::vectorSize() const {
  return propagatedSeed() ? propagatedSeed()->state().vectorSize() : 0; 
}

SeedUpdater *
SeedUpdater::ManagerImpl::createNewInstance(const String & key) {
  return new SeedUpdater(String("Update Seed ") + key);
}

} /* end namespace Pita */
