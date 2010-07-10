#include "OrthoDynamStateBasis.h"

namespace Pita {

OrthoDynamStateBasis::OrthoDynamStateBasis(size_t vectorSize, const DynamOps * metric, double tol) :
  vectorSize_(vectorSize),
  metric_(metric),
  orthoBasis_(DynamStatePlainBasis::New(vectorSize)),
  dualBasis_(DynamStatePlainBasis::New(vectorSize)),
  tolerance_(tol)
{}

/* Interface */

void
OrthoDynamStateBasis::originalBasisIs(const DynamStateBasis & ob) {
  checkSize(ob);

  // Reset results
  orthoBasis_->stateBasisDel();
  dualBasis_->stateBasisDel();
 
  addBasis(ob);
}

void
OrthoDynamStateBasis::originalBasisInc(const DynamStateBasis & ob) {
  checkSize(ob); 
  addBasis(ob);
}

/* Implementation */

void
OrthoDynamStateBasis::checkSize(const DynamStateBasis & basis) {
  if (basis.vectorSize() != vectorSize())
    throw RangeException();
}

void
OrthoDynamStateBasis::addBasis(const DynamStateBasis & originalBasis) {
  // Modified Gram-Schmidt procedure
  for (DynamStateBasis::IteratorConst it = originalBasis.state(); it; ++it) {
    DynamState currentState = *it;

    // Substract projections 
    size_t statesInBasis = dualBasis_->stateCount();
    for (size_t j = 0; j < statesInBasis; ++j) {
      double coef = currentState * dualBasis_->state(j);
      currentState.linAdd(-coef, orthoBasis_->state(j));
    }

    DynamState dualState = mult(metric(), currentState);

    // Skip if linearly dependent
    double current_norm = std::sqrt(currentState * dualState);
    if (current_norm < tolerance())
      continue;

    // Normalize and add
    double inverse_current_norm = 1.0 / current_norm;
    currentState.displacement() *= inverse_current_norm;
    currentState.velocity() *= inverse_current_norm;
    orthoBasis_->lastStateIs(currentState);
    dualState.displacement() *= inverse_current_norm;
    dualState.velocity() *= inverse_current_norm;
    dualBasis_->lastStateIs(dualState);
  } 
}

} /* end namespace Pita */
