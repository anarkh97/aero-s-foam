#ifndef ROM_GALERKINPROJECTIONSOLVER_H
#define ROM_GALERKINPROJECTIONSOLVER_H

#include "PodProjectionSolver.h"

#include "CholeskyUtils.h"

#include <Math.d/FullSquareMatrix.h>

#include <cstddef>
#include <cassert>

namespace Rom {

template <typename Scalar>
class GenGalerkinProjectionSolver : public GenPodProjectionSolver<Scalar> {
public:
  GenGalerkinProjectionSolver(Connectivity *cn, DofSetArray *dsa, ConstrainedDSA *c_dsa);

private:
  GenFullSquareMatrix<Scalar> reducedMatrix_;
  bool rhsIsprojected_;
 
  // Overriden 
  virtual void resetSolver(int vCount, int vSize);
  virtual void assembleAndFactorReducedSystem();
  virtual void projectRhs(const GenVector<Scalar> &);
  virtual double getReducedRhsNorm() const;
  virtual void solveReducedSystem(GenVector<Scalar> &);
 
  // Implementation 
  void performFactor();
  void performSolve();
};

template <typename Scalar>
GenGalerkinProjectionSolver<Scalar>::GenGalerkinProjectionSolver(Connectivity *cn,
                                                                 DofSetArray *dsa,
                                                                 ConstrainedDSA *c_dsa):
  GenPodProjectionSolver<Scalar>(cn, dsa, c_dsa),
  reducedMatrix_(),
  rhsIsprojected_(false)
{}

template <typename Scalar>
void
GenGalerkinProjectionSolver<Scalar>::resetSolver(int vCount, int) {
  reducedMatrix_.setSize(vCount);
  rhsIsprojected_ = false;
}

template <typename Scalar>
void
GenGalerkinProjectionSolver<Scalar>::assembleAndFactorReducedSystem() {
  for (int row = 0; row < this->basisSize(); ++row) {
    const GenVector<Scalar> &action = this->lastReducedMatrixAction()[row];
    for (int col = row; col < this->basisSize(); ++col) {
      reducedMatrix_[row][col] = action * this->projectionBasis()[col];
    }
  }

  performFactor();
  rhsIsprojected_ = false;
}

template <typename Scalar>
void
GenGalerkinProjectionSolver<Scalar>::projectRhs(const GenVector<Scalar> &rhs) {
  reduce(this->projectionBasis(), rhs, this->getReducedSolution());
  rhsIsprojected_ = true;
}

template <typename Scalar>
double
GenGalerkinProjectionSolver<Scalar>::getReducedRhsNorm() const {
  return this->lastReducedSolution().norm();
}

template <typename Scalar>
void
GenGalerkinProjectionSolver<Scalar>::solveReducedSystem(GenVector<Scalar> &rhs) {
  if (!rhsIsprojected_) {
    reduce(this->projectionBasis(), rhs, this->getReducedSolution());
  }
  
  performSolve();
  
  rhsIsprojected_ = false;
}

template <typename Scalar>
void
GenGalerkinProjectionSolver<Scalar>::performFactor() {
  cholesky_factor_upper(reducedMatrix_);
}

template <typename Scalar>
void
GenGalerkinProjectionSolver<Scalar>::performSolve() {
  cholesky_solve_upper(reducedMatrix_, this->getReducedSolution().data());
}

typedef GenGalerkinProjectionSolver<double> GalerkinProjectionSolver;

} /* end namespace Rom */

#endif /* ROM_GALERKINPROJECTIONSOLVER_H */
