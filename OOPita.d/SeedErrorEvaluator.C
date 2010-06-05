#include "SeedErrorEvaluator.h"

#include "DynamStateOps.h"

namespace Pita {

SeedErrorEvaluator::SeedErrorEvaluator(const DynamOps * metric, const Seed * seedError, const Seed * referenceSeed) :
  SharedState<DynamState>::NotifieeConst(seedError),
  metric_(metric),
  referenceSeed_(referenceSeed)
{}

void
SeedErrorEvaluator::onIteration() {
  double jumpEnergy = energy(metric(), seedError()->state());
  double actualSeedEnergy = energy(metric(), referenceSeed()->state());
  double energyRatio = (actualSeedEnergy != 0.0) ? jumpEnergy / actualSeedEnergy : -1.0;

  // HACK
  log() << "*** " << seedError()->name() << "_" << seedError()->iteration() << " relative energy = jumpEnergy / actualSeedEnergy = " << energyRatio << "\n";
}

SeedErrorEvaluator::Manager::Manager(const DynamOps * defaultMetric) :
  defaultMetric_(defaultMetric)
{}

SeedErrorEvaluator *
SeedErrorEvaluator::Manager::createNewInstance(const SeedErrorEvaluator::Manager::KeyType & key) {
  return new SeedErrorEvaluator(defaultMetric(), key, NULL);
}

} /* end namespace Pita */
