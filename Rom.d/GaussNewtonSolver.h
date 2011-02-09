#ifndef ROM_GAUSSNETONSOLVER_H
#define ROM_GAUSSNETONSOLVER_H

#include <Solvers.d/Solver.h>
#include <Math.d/DBSparseMatrix.h>

#include "VecBasis.h"
#include "LeastSquaresSolver.h"

#include "BasisOps.h"

#include <algorithm>
#include <cstddef>
#include <stdexcept>

template <typename Scalar>
class GenGaussNewtonSolver : public GenSolver<Scalar>, public GenDBSparseMatrix<Scalar> {
public:
  GenGaussNewtonSolver(Connectivity *cn, DofSetArray *dsa, ConstrainedDSA *c_dsa);

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

private:
  int basisSize_;
  const GenVecBasis<Scalar> *projectionBasis_;
  GenLeastSquaresSolver<Scalar> lsSolver_;
  
  GenVecBasis<Scalar> matrixAction_;
  GenVector<Scalar> reducedSolution_;

  void fillRhsBuffer(const GenVector<Scalar> &);
  void validateRhs(const GenVector<Scalar> &);

  // Disallow copy and assignment
  GenGaussNewtonSolver(const GenGaussNewtonSolver<Scalar> &);
  GenGaussNewtonSolver<Scalar> &operator=(const GenGaussNewtonSolver<Scalar> &);
};

template <typename Scalar>
GenGaussNewtonSolver<Scalar>::GenGaussNewtonSolver(Connectivity *cn,
                                                   DofSetArray *dsa,
                                                   ConstrainedDSA *c_dsa):
  GenDBSparseMatrix<Scalar>(cn, dsa, c_dsa),
  basisSize_(0),
  projectionBasis_(NULL),
  reducedSolution_(0)
{}

template <typename Scalar>
long
GenGaussNewtonSolver<Scalar>::size() {
  return GenDBSparseMatrix<Scalar>::size();
}

template <typename Scalar>
int
GenGaussNewtonSolver<Scalar>::neqs() {
  return GenDBSparseMatrix<Scalar>::neqs();
}

template <typename Scalar>
void
GenGaussNewtonSolver<Scalar>::zeroAll() {
  GenDBSparseMatrix<Scalar>::zeroAll();
}

template <typename Scalar>
void
GenGaussNewtonSolver<Scalar>::add(GenFullSquareMatrix<Scalar> &elMat, int *dofs) {
  GenDBSparseMatrix<Scalar>::add(elMat, dofs);
}

template <typename Scalar>
void
GenGaussNewtonSolver<Scalar>::addDiscreteMass(int dof, Scalar dmass) {
  GenDBSparseMatrix<Scalar>::addDiscreteMass(dof, dmass);
}

template <typename Scalar>
void
GenGaussNewtonSolver<Scalar>::projectionBasisIs(const GenVecBasis<Scalar> &reducedBasis) {
  if (reducedBasis.vectorSize() != neqs()) {
    throw std::domain_error("Vectors of the reduced basis have the wrong size");
  }
  
  lsSolver_.problemSizeIs(reducedBasis.vectorSize(), reducedBasis.vectorCount());
  
  GenVecBasis<Scalar> temp(reducedBasis.vectorCount(), reducedBasis.vectorSize());
  temp.swap(matrixAction_);
  reducedSolution_.reset(reducedBasis.vectorCount(), Scalar());

  projectionBasis_ = &reducedBasis;
  basisSize_ = reducedBasis.vectorCount();
}

template <typename Scalar>
void
GenGaussNewtonSolver<Scalar>::factor() {
  lsSolver_.statusIs(LeastSquares::READY);
  
  for (int col = 0; col < lsSolver_.unknownCount(); ++col) {
    GenDBSparseMatrix<Scalar>::mult((*projectionBasis_)[col], matrixAction_[col]);
  }
    
  std::copy(&matrixAction_[0][0],
            &matrixAction_[0][0] + lsSolver_.equationCount() * lsSolver_.unknownCount(),
            lsSolver_.matrixBuffer());

  lsSolver_.statusIs(LeastSquares::FACTORED);
}

template <typename Scalar>
inline
void
GenGaussNewtonSolver<Scalar>::fillRhsBuffer(const GenVector<Scalar> &rhs) {
  std::copy(rhs.data(), rhs.data() + rhs.size(), lsSolver_.rhsBuffer(0));
}

template <typename Scalar>
inline
void
GenGaussNewtonSolver<Scalar>::validateRhs(const GenVector<Scalar> &rhs) {
  if (rhs.size() != neqs()) {
    throw std::domain_error("Rhs has the wrong size");
  }
}

template <typename Scalar>
void
GenGaussNewtonSolver<Scalar>::reSolve(GenVector<Scalar> &rhs) {
  validateRhs(rhs);

  assert(lsSolver_.status() == LeastSquares::FACTORED ||
         lsSolver_.status() == LeastSquares::PROJECTED);
  
  if (lsSolver_.status() == LeastSquares::FACTORED) {
    fillRhsBuffer(rhs);
  }

  lsSolver_.statusIs(LeastSquares::SOLVED);
  std::copy(lsSolver_.rhsBuffer(0),
            lsSolver_.rhsBuffer(0) + lsSolver_.unknownCount(),
            reducedSolution_.data());
  lsSolver_.statusIs(LeastSquares::FACTORED);
  
  expand(*projectionBasis_, reducedSolution_, rhs);
}

template <typename Scalar>
double
GenGaussNewtonSolver<Scalar>::projectAndComputeNorm(const GenVector<Scalar> &rhs) {
  validateRhs(rhs);

  assert(lsSolver_.status() == LeastSquares::FACTORED);

  fillRhsBuffer(rhs);
  lsSolver_.statusIs(LeastSquares::PROJECTED);

  const GenVector<Scalar> components(lsSolver_.rhsBuffer(0), lsSolver_.unknownCount(), false);
  return components.norm();
}

typedef GenGaussNewtonSolver<double> GaussNewtonSolver;

#endif /* ROM_GAUSSNETONSOLVER_H */
