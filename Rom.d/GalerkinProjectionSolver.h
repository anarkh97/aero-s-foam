#ifndef ROM_GALERKINPROJECTIONSOLVER_H
#define ROM_GALERKINPROJECTIONSOLVER_H

#include <Solvers.d/Solver.h>
#include <Math.d/DBSparseMatrix.h>

#include <Math.d/VectorSet.h>

#include <Math.d/FullSquareMatrix.h>

#include <memory>
#include <cstddef>

#include <cassert>

template <typename Scalar>
class GenGalerkinProjectionSolver : public GenSolver<Scalar>, public GenDBSparseMatrix<Scalar> {
public:
  GenGalerkinProjectionSolver(Connectivity *cn, DofSetArray *dsa, ConstrainedDSA *c_dsa);

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

  // Projection basis
  int basisSize() const { return basisSize_; }
  const GenVectorSet<Scalar> &projectionBasis() const { return *projectionBasis_; }
  void projectionBasisIs(const GenVectorSet<Scalar> &); // Passed object must be kept alive by owner

  // Data collection
  const GenVector<Scalar> &lastReducedSolution() const { return reducedRhs_; }
  const GenVectorSet<Scalar> &lastReducedMatrixAction() const { return *matrixAction_; }

private:
  void factorReducedMatrix();
  void solveReducedRhs();

  int basisSize_;
  const GenVectorSet<Scalar> *projectionBasis_;
  
  GenFullSquareMatrix<Scalar> reducedMatrix_;

  std::auto_ptr<GenVectorSet<Scalar> > matrixAction_;
  GenVector<Scalar> reducedRhs_;

  // Disallow copy and assignment
  GenGalerkinProjectionSolver(const GenGalerkinProjectionSolver<Scalar> &);
  GenGalerkinProjectionSolver<Scalar> &operator=(const GenGalerkinProjectionSolver<Scalar> &);
};

template <typename Scalar>
GenGalerkinProjectionSolver<Scalar>::GenGalerkinProjectionSolver(Connectivity *cn,
                                                                 DofSetArray *dsa,
                                                                 ConstrainedDSA *c_dsa):
  GenDBSparseMatrix<Scalar>(cn, dsa, c_dsa),
  basisSize_(0),
  projectionBasis_(NULL),
  reducedMatrix_(),
  matrixAction_(new GenVectorSet<Scalar>()),
  reducedRhs_(0)
{
  projectionBasis_ = matrixAction_.get();
}

template <typename Scalar>
long
GenGalerkinProjectionSolver<Scalar>::size() {
  return GenDBSparseMatrix<Scalar>::size();
}

template <typename Scalar>
int
GenGalerkinProjectionSolver<Scalar>::neqs() {
  return GenDBSparseMatrix<Scalar>::neqs();
}

template <typename Scalar>
void
GenGalerkinProjectionSolver<Scalar>::zeroAll() {
  GenDBSparseMatrix<Scalar>::zeroAll();
}

template <typename Scalar>
void
GenGalerkinProjectionSolver<Scalar>::add(GenFullSquareMatrix<Scalar> &elMat, int *dofs) {
  GenDBSparseMatrix<Scalar>::add(elMat, dofs);
}

template <typename Scalar>
void
GenGalerkinProjectionSolver<Scalar>::addDiscreteMass(int dof, Scalar dmass) {
  GenDBSparseMatrix<Scalar>::addDiscreteMass(dof, dmass);
}

template <typename Scalar>
void
GenGalerkinProjectionSolver<Scalar>::projectionBasisIs(const GenVectorSet<Scalar> &basis) {
  assert(neqs() == basis.size()); // TODO: Exception

  const int newBasisSize = basis.numVec();
  std::auto_ptr<GenVectorSet<Scalar> > newMatrixAction(new GenVectorSet<Scalar>(newBasisSize, neqs()));
  GenVector<Scalar> newReducedRhs(newBasisSize);

  reducedMatrix_.setSize(newBasisSize);

  basisSize_ = newBasisSize;
  projectionBasis_ = &basis;
  matrixAction_ = newMatrixAction;
  reducedRhs_.swap(newReducedRhs);
}

template <typename Scalar>
void
GenGalerkinProjectionSolver<Scalar>::factor() {
  for (int row = 0; row < basisSize(); ++row) {
    GenVector<Scalar> &action = (*matrixAction_)[row];
    GenDBSparseMatrix<Scalar>::mult((*projectionBasis_)[row], action);
    for (int col = row; col < basisSize(); ++col) {
      reducedMatrix_[row][col] = action * (*projectionBasis_)[col];
    }
  }

  factorReducedMatrix();
}

template <typename Scalar>
void
GenGalerkinProjectionSolver<Scalar>::reSolve(GenVector<Scalar> &rhs) {
  assert(neqs() == rhs.size()); // TODO: Exception
  
  assert(basisSize() == reducedMatrix_.dim());
  assert(basisSize() == reducedRhs_.size());
  
  for (int i = 0; i < basisSize(); ++i) {
    reducedRhs_[i] = (*projectionBasis_)[i] * rhs;
  }

  solveReducedRhs();

  rhs.zero();
  for (int i = 0; i < basisSize(); ++i) {
    rhs.linAdd(reducedRhs_[i], (*projectionBasis_)[i]);
  }
}

typedef GenGalerkinProjectionSolver<double> GalerkinProjectionSolver;

#endif /* ROM_GALERKINPROJECTIONSOLVER_H */
