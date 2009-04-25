#include "AlternatingIntegratorPropagator.h"

namespace Pita {

AlternatingIntegratorPropagator::AlternatingIntegratorPropagator(DynamTimeIntegrator * integrator) :
  IntegratorPropagator(integrator),
  secondaryPropagator_(NULL)
{}

void
AlternatingIntegratorPropagator::initialStateIs(const DynamState & initialState) {
  if (secondaryPropagator_) {
    setInitialState(initialState);
    initialStateNotify();

    if (integrator()) {
      integrator()->currentSliceIs(this->sliceRank());
      integrator()->initialConditionIs(initialState, initialTime());
      unsigned int steps = timeStepCount().value();
      for (unsigned int s = 0; s < steps; ++s) {
        integrator()->timeStepCountInc();
        secondaryPropagator_->timeStepCountInc();
      }
      setFinalState(integrator()->currentState());
    } else {
      setFinalState(initialState);
    }

    finalStateNotify();
  } else {
    // Falls back to unspecialized propagation
    IntegratorPropagator::initialStateIs(initialState);
  }
}

} // end namespace Pita
