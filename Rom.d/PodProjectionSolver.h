#ifndef ROM_PODPROJECTIONSOLVER_H
#define ROM_PODPROJECTIONSOLVER_H

#include <Solvers.d/Solver.h>
#include <Math.d/DBSparseMatrix.h>

#include "VecBasis.h"
#include "BasisOps.h"

#include <cstddef>
#include <stdexcept>

template <typename Scalar>
class GenPodProjectionSolver : public GenSolver<Scalar>, public GenDBSparseMatrix<Scalar> {
public:
  GenPodProjectionSolver(Connectivity *cn, DofSetArray *dsa, ConstrainedDSA *c_dsa);

  // Pure virtual function implementations
  virtual long size();
  virtual int neqs();

  // Full-order matrix assembly
  virtual void zeroAll();
  virtual void add(GenFullSquareMatrix<Scalar> &, int *);
  virtual void addDiscreteMass(int, Scalar);

  // Solution
  virtual void factor();
  virtual void reSolve(GenVector<Scalar> &rhs);
  double projectAndComputeNorm(const GenVector<Scalar> &rhs); // next reSolve must use same rhs

  // Reduced basis parameters
  int basisSize() const { return basisSize_; }
  const GenVecBasis<Scalar> &projectionBasis() const { return *projectionBasis_; }
  void projectionBasisIs(const GenVecBasis<Scalar> &); // Passed objects must be kept alive by owner
  
  // Data collection
  const GenVector<Scalar> &lastReducedSolution() const { return reducedSolution_; }
  const GenVecBasis<Scalar> &lastReducedMatrixAction() const { return matrixAction_; }

protected:
  GenVector<Scalar> &getReducedSolution() { return reducedSolution_; }

private:
  int basisSize_;
  const GenVecBasis<Scalar> *projectionBasis_;
  
  GenVecBasis<Scalar> matrixAction_;
  GenVector<Scalar> reducedSolution_;
  
  void validateRhs(const GenVector<Scalar> &);

  virtual void resetSolver(int vCount, int vSize) = 0;
  virtual void assembleAndFactorReducedSystem() = 0;
  virtual void projectRhs(const GenVector<Scalar> &) = 0;
  virtual double getReducedRhsNorm() const = 0;
  virtual void solveReducedSystem(GenVector<Scalar> &) = 0;

  // Disallow copy and assignment
  GenPodProjectionSolver(const GenPodProjectionSolver<Scalar> &);
  GenPodProjectionSolver<Scalar> &operator=(const GenPodProjectionSolver<Scalar> &);
};

template <typename Scalar>
GenPodProjectionSolver<Scalar>::GenPodProjectionSolver(Connectivity *cn,
                                                   DofSetArray *dsa,
                                                   ConstrainedDSA *c_dsa):
  GenDBSparseMatrix<Scalar>(cn, dsa, c_dsa),
  basisSize_(0),
  projectionBasis_(NULL),
  matrixAction_(0, 0),
  reducedSolution_(0)
{
  projectionBasis_ = &matrixAction_;
}

template <typename Scalar>
long
GenPodProjectionSolver<Scalar>::size() {
  return GenDBSparseMatrix<Scalar>::size();
}

template <typename Scalar>
int
GenPodProjectionSolver<Scalar>::neqs() {
  return GenDBSparseMatrix<Scalar>::neqs();
}

template <typename Scalar>
void
GenPodProjectionSolver<Scalar>::zeroAll() {
  GenDBSparseMatrix<Scalar>::zeroAll();
}

template <typename Scalar>
void
GenPodProjectionSolver<Scalar>::add(GenFullSquareMatrix<Scalar> &elMat, int *dofs) {
  GenDBSparseMatrix<Scalar>::add(elMat, dofs);
}

template <typename Scalar>
void
GenPodProjectionSolver<Scalar>::addDiscreteMass(int dof, Scalar dmass) {
  GenDBSparseMatrix<Scalar>::addDiscreteMass(dof, dmass);
}

template <typename Scalar>
void
GenPodProjectionSolver<Scalar>::projectionBasisIs(const GenVecBasis<Scalar> &reducedBasis) {
  if (reducedBasis.vectorSize() != neqs()) {
    throw std::domain_error("Vectors of the reduced basis have the wrong size");
  }

  GenVecBasis<Scalar> newMatrixAction(reducedBasis.vectorCount(), reducedBasis.vectorSize());
  GenVector<Scalar> newReducedSolution(reducedBasis.vectorCount());

  resetSolver(reducedBasis.vectorCount(), reducedBasis.vectorSize());
  
  swap(matrixAction_, newMatrixAction);
  reducedSolution_.swap(newReducedSolution);
  projectionBasis_ = &reducedBasis;
  basisSize_ = reducedBasis.vectorCount();
}

template <typename Scalar>
inline
void
GenPodProjectionSolver<Scalar>::validateRhs(const GenVector<Scalar> &rhs) {
  if (rhs.size() != neqs()) {
    throw std::domain_error("Rhs has the wrong size");
  }
}

template <typename Scalar>
void
GenPodProjectionSolver<Scalar>::factor() {
  for (int col = 0; col < basisSize_; ++col) {
    GenDBSparseMatrix<Scalar>::mult(projectionBasis()[col], matrixAction_[col]);
  }

  assembleAndFactorReducedSystem();
}
   
template <typename Scalar>
void
GenPodProjectionSolver<Scalar>::reSolve(GenVector<Scalar> &rhs) {
  validateRhs(rhs);
  solveReducedSystem(rhs);
  expand(projectionBasis(), lastReducedSolution(), rhs);
}

template <typename Scalar>
double
GenPodProjectionSolver<Scalar>::projectAndComputeNorm(const GenVector<Scalar> &rhs) {
  validateRhs(rhs);
  projectRhs(rhs);
  return getReducedRhsNorm();
}

typedef GenPodProjectionSolver<double> PodProjectionSolver;

#endif /* ROM_PODPROJECTIONSOLVER_H */
