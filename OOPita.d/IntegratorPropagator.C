#include "IntegratorPropagator.h"
#include "DynamTimeIntegrator.h"

namespace Pita {
  
IntegratorPropagator::IntegratorPropagator(DynamTimeIntegrator * integrator) :
  DynamPropagator(integrator ? integrator->vectorSize() : 0),
  integrator_(integrator),
  initialTime_(Seconds(0.0)),
  timeStepCount_(TimeStepCount(1))
{}

void
IntegratorPropagator::initialStateIs(const DynamState & initialState) {
  setInitialState(initialState);
  initialStateNotify();

  if (integrator()) {
    integrator()->currentSliceIs(this->sliceRank());
    integrator()->initialConditionIs(initialState, initialTime());
    integrator()->timeStepCountInc(this->timeStepCount());
    setFinalState(integrator()->currentState());
  } else {
    setFinalState(initialState);
  }

  finalStateNotify();
}

} // end namespace Pita
