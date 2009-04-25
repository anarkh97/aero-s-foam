#include "LinearizedPropagator.h"
#include "NlDynamTimeIntegrator.h" 
#include "PitaNonLinDynam.h"

namespace Pita {

LinearizedPropagator::LinearizedPropagator(const NlDynamTimeIntegrator * integrator) :
  DynamPropagator(integrator->vectorSize()),
  integrator_(integrator),
  probDesc_(dynamic_cast<const PitaNonLinDynamic*>(integrator->probDesc())),
  temp_(integrator->vectorSize()),
  midTimeDisp_(integrator->vectorSize())
{}

void
LinearizedPropagator::initialStateIs(const DynamState & initialState) {
  double delta = 0.5 * integrator_->timeStepSize().value();
  temp_ = initialState.displacement();
  temp_.linAdd(delta, initialState.velocity());
  probDesc_->zeroRotDofs(temp_);
  const_cast<GenSparseMatrix<double>*>(const_cast<PitaNonLinDynamic*>(probDesc_)->getMassMatrix())->mult(temp_, midTimeDisp_);
  const_cast<PitaNonLinDynamic*>(probDesc_)->getSolver()->reSolve(midTimeDisp_);
  temp_.linC(-1.0, initialState.displacement(), 2.0, midTimeDisp_);
  DynamState finalState(initialState.vectorSize());
  finalState.displacement() = temp_;
  temp_ -= initialState.displacement();
  temp_ *= (1.0 / delta);
  temp_ -= initialState.velocity();
  probDesc_->zeroRotDofs(temp_);
  finalState.velocity() = temp_;

  // Commit changes
  setInitialState(initialState);
  initialStateNotify();
  setFinalState(finalState);
  finalStateNotify();
}

} // end namespace Pita
