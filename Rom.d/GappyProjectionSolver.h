#ifndef ROM_GAPPYPROJECTIONSOLVER_H
#define ROM_GAPPYPROJECTIONSOLVER_H

#include <Solvers.d/Solver.h>
#include <Math.d/DBSparseMatrix.h>

#include "VecBasis.h"
#include "NodalRestrictionMapping.h"
#include "LeastSquaresSolver.h"

#include "BasisOps.h"

#include <cstddef>
#include <stdexcept>
#include <cassert>

namespace Rom {

template <typename Scalar>
class GenGappyProjectionSolver : public GenSolver<Scalar>, public GenDBSparseMatrix<Scalar> {
public:
  GenGappyProjectionSolver(Connectivity *cn, DofSetArray *dsa, ConstrainedDSA *c_dsa);

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
  int reducedBasisSize() const { return reducedBasisSize_; }       // n_y
  int projectionBasisSize() const { return projectionBasisSize_; } // n_J
  int vectorSize() const { return neqs(); }                        // n_j (mesh size)
  int sampleSize() const { return sampleSize_; }                   // n_i (sample dofs)

  // Passed objects must be kept alive by owner
  void systemApproximationIs(const NodalRestrictionMapping &sampleMapping,   // n_i <-> n_j 
                             const GenVecBasis<Scalar> &reducedBasis,        // Phi_y: n_y * n_j
                             const GenVecBasis<Scalar> &jacobianProjection,  // A:     n_J * n_i
                             const GenVecBasis<Scalar> &residualProjection); // B:     n_J * n_i

private:
  int reducedBasisSize_;
  int projectionBasisSize_;
  int sampleSize_;

  const NodalRestrictionMapping *sampleMapping_;

  const GenVecBasis<Scalar> *reducedBasis_;
  const GenVecBasis<Scalar> *jacobianProjection_;
  const GenVecBasis<Scalar> *residualProjection_;

  GenLeastSquaresSolver<Scalar> lsSolver_;
  
  void fillRhsBuffer(const GenVector<Scalar> &);
  void validateRhs(const GenVector<Scalar> &);
  
  // Disallow copy and assignment
  GenGappyProjectionSolver(const GenGappyProjectionSolver<Scalar> &);
  GenGappyProjectionSolver<Scalar> &operator=(const GenGappyProjectionSolver<Scalar> &);
};

template <typename Scalar>
GenGappyProjectionSolver<Scalar>::GenGappyProjectionSolver(Connectivity *cn,
                                                           DofSetArray *dsa,
                                                           ConstrainedDSA *c_dsa):
  GenDBSparseMatrix<Scalar>(cn, dsa, c_dsa),
  reducedBasisSize_(0),
  projectionBasisSize_(0),
  sampleSize_(0),
  sampleMapping_(NULL),
  reducedBasis_(NULL),
  jacobianProjection_(NULL),
  residualProjection_(NULL)
{}

template <typename Scalar>
long
GenGappyProjectionSolver<Scalar>::size() {
  return GenDBSparseMatrix<Scalar>::size();
}

template <typename Scalar>
int
GenGappyProjectionSolver<Scalar>::neqs() {
  return GenDBSparseMatrix<Scalar>::neqs();
}

template <typename Scalar>
void
GenGappyProjectionSolver<Scalar>::zeroAll() {
  GenDBSparseMatrix<Scalar>::zeroAll();
}

template <typename Scalar>
void
GenGappyProjectionSolver<Scalar>::add(GenFullSquareMatrix<Scalar> &elMat, int *dofs) {
  GenDBSparseMatrix<Scalar>::add(elMat, dofs);
}

template <typename Scalar>
void
GenGappyProjectionSolver<Scalar>::addDiscreteMass(int dof, Scalar dmass) {
  GenDBSparseMatrix<Scalar>::addDiscreteMass(dof, dmass);
}

template <typename Scalar>
void
GenGappyProjectionSolver<Scalar>::systemApproximationIs(const NodalRestrictionMapping &sampleMapping,
                                                        const GenVecBasis<Scalar> &reducedBasis,
                                                        const GenVecBasis<Scalar> &jacobianProjection,
                                                        const GenVecBasis<Scalar> &residualProjection) {
  assert(sampleMapping.originInfo() == reducedBasis.vectorSize());
  assert(sampleMapping.restrictedInfo() == jacobianProjection.vectorSize());
  assert(sampleMapping.restrictedInfo() == residualProjection.vectorSize());
  assert(jacobianProjection.vectorCount() == residualProjection.vectorCount());

  assert(reducedBasis.vectorCount() <= jacobianProjection.vectorCount());
  assert(jacobianProjection.vectorCount() <= sampleMapping.restrictedInfo());

  reducedBasisSize_ = reducedBasis.vectorCount();
  projectionBasisSize_ = jacobianProjection.vectorCount();
  sampleSize_ = sampleMapping.restrictedInfo();
 
  lsSolver_.problemSizeIs(projectionBasisSize_, reducedBasisSize_);

  sampleMapping_ = &sampleMapping;
  reducedBasis_  = &reducedBasis;
  jacobianProjection_ = &jacobianProjection;
  residualProjection_ = &residualProjection;
}

template <typename Scalar>
void
GenGappyProjectionSolver<Scalar>::factor() {
  lsSolver_.statusIs(LeastSquares::READY);
  
  GenVector<Scalar> matrixAction(sampleMapping_->originInfo());
  for (int col = 0; col < reducedBasisSize(); ++col) {
    GenDBSparseMatrix<Scalar>::mult((*reducedBasis_)[col], matrixAction);
    for (int row = 0; row < projectionBasisSize(); ++row) {
      lsSolver_.matrixEntry(row, col) = sampleMapping_->dotProduct(matrixAction, (*jacobianProjection_)[row]);
    }
  }

  lsSolver_.statusIs(LeastSquares::FACTORED);
}

template <typename Scalar>
inline
void
GenGappyProjectionSolver<Scalar>::fillRhsBuffer(const GenVector<Scalar> &rhs) {
  for (int row = 0; row < projectionBasisSize(); ++row) {
    lsSolver_.rhsEntry(row) = sampleMapping_->dotProduct(rhs, (*residualProjection_)[row]);
  }
}

template <typename Scalar>
inline
void
GenGappyProjectionSolver<Scalar>::validateRhs(const GenVector<Scalar> &rhs) {
  if (rhs.size() != neqs()) {
    throw std::domain_error("Rhs has the wrong size");
  }
}

template <typename Scalar>
void
GenGappyProjectionSolver<Scalar>::reSolve(GenVector<Scalar> &rhs) {
  validateRhs(rhs);

  assert(lsSolver_.status() == LeastSquares::FACTORED ||
         lsSolver_.status() == LeastSquares::PROJECTED);
  
  if (lsSolver_.status() == LeastSquares::FACTORED) {
    fillRhsBuffer(rhs);
  }

  lsSolver_.statusIs(LeastSquares::SOLVED);
  GenVector<Scalar> reducedSolution(lsSolver_.rhsBuffer(0), lsSolver_.unknownCount(), false);
  lsSolver_.statusIs(LeastSquares::FACTORED);
  
  expand(*reducedBasis_, reducedSolution, rhs);
}

template <typename Scalar>
double
GenGappyProjectionSolver<Scalar>::projectAndComputeNorm(const GenVector<Scalar> &rhs) {
  validateRhs(rhs);

  assert(lsSolver_.status() == LeastSquares::FACTORED ||
         lsSolver_.status() == LeastSquares::PROJECTED);

  lsSolver_.statusIs(LeastSquares::FACTORED);
  fillRhsBuffer(rhs);
  lsSolver_.statusIs(LeastSquares::PROJECTED);

  const GenVector<Scalar> components(lsSolver_.rhsBuffer(0), lsSolver_.unknownCount(), false);
  return components.norm();
}

typedef GenGappyProjectionSolver<double> GappyProjectionSolver;

} /* end namespace Rom */

#endif /* ROM_GAPPYPROJECTIONSOLVER_H */
