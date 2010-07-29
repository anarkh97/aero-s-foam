#include "SeedDifferenceEvaluator.h"

#include "DynamStateOps.h"
#include "NlDynamOps.h"
#include <Math.d/SparseMatrix.h>

namespace Pita {

void
SeedDifferenceEvaluator::onIteration() {
  double refMagnitude = getRefMagnitude();
  double magnitudeRatio = (refMagnitude != 0.0) ? (getDiffMagnitude() / refMagnitude) : -1.0;
  
  // HACK, find a better outlet
  log() << "JUMP " << seedDifference()->name() << " at iter " << seedDifference()->iteration()
        << " -- Relative error = " << magnitudeRatio << "\n";
}

LinSeedDifferenceEvaluator::LinSeedDifferenceEvaluator(const DynamOps * metric, const Seed * seedDifference, const Seed * referenceSeed) :
  SeedDifferenceEvaluator(seedDifference, referenceSeed),
  metric_(metric)
{}

double
LinSeedDifferenceEvaluator::getRefMagnitude() const {
  return energy(metric(), referenceSeed()->state());
}

double
LinSeedDifferenceEvaluator::getDiffMagnitude() const {
  return energy(metric(), seedDifference()->state());
}

LinSeedDifferenceEvaluator::Manager::Manager(const DynamOps * defaultMetric) :
  defaultMetric_(defaultMetric)
{}

LinSeedDifferenceEvaluator::Ptr
LinSeedDifferenceEvaluator::Manager::createNewInstance(const LinSeedDifferenceEvaluator::Manager::KeyType & key) {
  return new LinSeedDifferenceEvaluator(defaultMetric(), key, NULL);
}

NonLinSeedDifferenceEvaluator::NonLinSeedDifferenceEvaluator(PitaNonLinDynamic * probDesc,
                                                             const Seed * seedDifference,
                                                             const Seed * referenceSeed,
                                                             Seconds refTime) :
  SeedDifferenceEvaluator(seedDifference, referenceSeed),
  probDesc_(probDesc),
  refTime_(refTime) 
{}

double
NonLinSeedDifferenceEvaluator::getRefMagnitude() const {
  Vector dualVelocity(referenceSeed()->state().vectorSize());
  SparseMatrix * massMatrix = const_cast<SparseMatrix *>(const_cast<PitaNonLinDynamic *>(probDesc_)->getMassMatrix());
  massMatrix->mult(referenceSeed()->state().velocity(), dualVelocity);

  double refKineticEnergy = 0.5 * (dualVelocity * referenceSeed()->state().velocity());
  double refInternalEnergy = probDesc_->internalEnergy(referenceSeed()->state().displacement());
  return refKineticEnergy + refInternalEnergy;
}

double
NonLinSeedDifferenceEvaluator::getDiffMagnitude() const {
  NlDynamOps::Ptr metric = NlDynamOps::New(probDesc_);
  metric->displacementIs(referenceSeed()->state().displacement(), refTime_);
  return energy(metric.ptr(), seedDifference()->state());
}

NonLinSeedDifferenceEvaluator::Ptr
NonLinSeedDifferenceEvaluator::Manager::createNewInstance(const NonLinSeedDifferenceEvaluator::Manager::KeyType & key) {
  return new NonLinSeedDifferenceEvaluator(probDesc_, key, NULL);
}

} /* end namespace Pita */
