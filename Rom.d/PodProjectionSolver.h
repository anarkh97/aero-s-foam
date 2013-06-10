#ifndef ROM_PODPROJECTIONSOLVER_H
#define ROM_PODPROJECTIONSOLVER_H

#include <Solvers.d/Solver.h>
#include <Math.d/DBSparseMatrix.h>

#include "VecBasis.h"
#include "BasisOps.h"

#include <cstddef>
#include <stdexcept>

namespace Rom {

template <typename Scalar>
class GenPodProjectionSolver {
public:
  // Solution
  virtual void factor() = 0;
  virtual void reSolve(GenVector<Scalar> &rhs) = 0;
  virtual double projectAndComputeNorm(const GenVector<Scalar> &rhs) = 0; // next reSolve must use same rhs

  // Reduced basis parameters
  virtual int basisSize() const = 0;
  virtual const GenVecBasis<Scalar> &projectionBasis() const = 0;
  virtual void projectionBasisIs(const GenVecBasis<Scalar> &) = 0; 
  virtual void projectionBasis2Is(const GenVecBasis<Scalar> &) {};

  // Data collection
  virtual const GenVector<Scalar> &lastReducedSolution() const = 0;
  virtual const GenVecBasis<Scalar> &lastReducedMatrixAction() const = 0;
};

typedef GenPodProjectionSolver<double> PodProjectionSolver;

template <typename Scalar>
class GenDBSparsePodProjectionSolver : public GenPodProjectionSolver<Scalar>, public GenSolver<Scalar>, public GenDBSparseMatrix<Scalar> {
public:
  GenDBSparsePodProjectionSolver(Connectivity *cn, DofSetArray *dsa, ConstrainedDSA *c_dsa);

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
  void projectionBasis2Is(const GenVecBasis<Scalar> &);
  
  // Data collection
  const GenVector<Scalar> &lastReducedSolution() const { return reducedSolution_; }
  const GenVecBasis<Scalar> &lastReducedMatrixAction() const { return matrixAction_; }

protected:
  GenVector<Scalar> &getReducedSolution() { return reducedSolution_; }

private:
  int basisSize_;
  const GenVecBasis<Scalar> *projectionBasis_;
  const GenVecBasis<Scalar> *projectionBasis2_;
  
  GenVecBasis<Scalar> matrixAction_;
  GenVector<Scalar> reducedSolution_;
  
  void validateRhs(const GenVector<Scalar> &);

  virtual void resetSolver(int vCount, int vSize) = 0;
  virtual void assembleAndFactorReducedSystem() = 0;
  virtual void projectRhs(const GenVector<Scalar> &) = 0;
  virtual double getReducedRhsNorm() const = 0;
  virtual void solveReducedSystem(GenVector<Scalar> &) = 0;

  // Disallow copy and assignment
  GenDBSparsePodProjectionSolver(const GenDBSparsePodProjectionSolver<Scalar> &);
  GenDBSparsePodProjectionSolver<Scalar> &operator=(const GenDBSparsePodProjectionSolver<Scalar> &);
};

template <typename Scalar>
GenDBSparsePodProjectionSolver<Scalar>::GenDBSparsePodProjectionSolver(Connectivity *cn,
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
GenDBSparsePodProjectionSolver<Scalar>::size() {
  return GenDBSparseMatrix<Scalar>::size();
}

template <typename Scalar>
int
GenDBSparsePodProjectionSolver<Scalar>::neqs() {
  return GenDBSparseMatrix<Scalar>::neqs();
}

template <typename Scalar>
void
GenDBSparsePodProjectionSolver<Scalar>::zeroAll() {
  GenDBSparseMatrix<Scalar>::zeroAll();
}

template <typename Scalar>
void
GenDBSparsePodProjectionSolver<Scalar>::add(GenFullSquareMatrix<Scalar> &elMat, int *dofs) {
  GenDBSparseMatrix<Scalar>::add(elMat, dofs);
}

template <typename Scalar>
void
GenDBSparsePodProjectionSolver<Scalar>::addDiscreteMass(int dof, Scalar dmass) {
  GenDBSparseMatrix<Scalar>::addDiscreteMass(dof, dmass);
}

template <typename Scalar>
void
GenDBSparsePodProjectionSolver<Scalar>::projectionBasisIs(const GenVecBasis<Scalar> &reducedBasis) {
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
void
GenDBSparsePodProjectionSolver<Scalar>::projectionBasis2Is(const GenVecBasis<Scalar> &reducedBasis2) {
  projectionBasis2_ = &reducedBasis2;
}


template <typename Scalar>
inline
void
GenDBSparsePodProjectionSolver<Scalar>::validateRhs(const GenVector<Scalar> &rhs) {
  if (rhs.size() != neqs()) {
    throw std::domain_error("Rhs has the wrong size");
  }
}

template <typename Scalar>
void
GenDBSparsePodProjectionSolver<Scalar>::factor() {
  for (int col = 0; col < basisSize_; ++col) {
    GenDBSparseMatrix<Scalar>::mult(projectionBasis()[col], matrixAction_[col]);
  }

  assembleAndFactorReducedSystem();
}
   
template <typename Scalar>
void
GenDBSparsePodProjectionSolver<Scalar>::reSolve(GenVector<Scalar> &rhs) {
/*
  validateRhs(rhs);
*/
  solveReducedSystem(rhs);
/*
  expand(projectionBasis(), lastReducedSolution(), rhs);
*/
}

template <typename Scalar>
double
GenDBSparsePodProjectionSolver<Scalar>::projectAndComputeNorm(const GenVector<Scalar> &rhs) {
/*
  validateRhs(rhs);
  projectRhs(rhs);
  return getReducedRhsNorm();
*/
  return rhs.norm();
}

typedef GenDBSparsePodProjectionSolver<double> DBSparsePodProjectionSolver;

} /* end namespace Rom */

#endif /* ROM_PODPROJECTIONSOLVER_H */
