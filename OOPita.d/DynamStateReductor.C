#include "DynamStateReductor.h"

namespace Pita {

DynamStateReductor::DynamStateReductor(const DynamStateBasis * reductionBasis,
                                       const RankDeficientSolver * solver) :
  reductionBasis_(reductionBasis),
  solver_(solver)
{}

Vector
DynamStateReductor::reducedComponents(const DynamState & is) const {
  assert(int(reducedBasisSize()) == solver()->matrixSize());

  if (is.vectorSize() != vectorSize()) {
    throw Fwk::RangeException("Dimension mismatch");
  }

  Vector result(reducedBasisSize(), 0.0);

  if (reducedBasisSize() != 0) {
    // Assemble relevant part of rhs
    for (int i = 0; i < solver_->factorRank(); ++i) {
      int index = solver()->factorPermutation(i);
      result[index] = is * reductionBasis()->state(index);
    }

    /*log() << "Rhs = ";
    for (int i = 0; i < reducedBasisSize(); ++i) {
      log() << result[i] << " ";
    }
    log() << "\n";*/
    
    // Perform in place resolution
    solver()->solution(result);

    /*log() << "Solution = ";
    for (int i = 0; i < reducedBasisSize(); ++i) {
      log() << result[i] << " ";
    }
    log() << "\n";*/
  }

  return result;
}

} /* end namespace Pita */
