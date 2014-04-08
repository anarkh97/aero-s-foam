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
  // Reduced-order matrix assembly
  virtual void addReducedMass(double Mcoef) = 0;

  // Solution
  virtual void factor() = 0;
  virtual void reSolve(GenVector<Scalar> &rhs) = 0;

  // Reduced basis parameters
  virtual int basisSize() const = 0;
  virtual const GenVecBasis<Scalar> &projectionBasis() const = 0;
  virtual void projectionBasisIs(const GenVecBasis<Scalar> &) = 0; 
  virtual void EmpiricalSolver() = 0; 
  virtual void addToReducedMatrix(const Eigen::Matrix<Scalar,Eigen::Dynamic,Eigen::Dynamic> &, double = 1.0) = 0;

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

  // Reduced matrix assembly
  void addReducedMass(double Mcoef) { Mcoef_ = Mcoef; }

  // Solution
  virtual void factor();
  virtual void reSolve(GenVector<Scalar> &rhs);

  // Reduced basis parameters
  int basisSize() const { return basisSize_; }
  const GenVecBasis<Scalar> &projectionBasis() const { return *projectionBasis_; }
  void projectionBasisIs(const GenVecBasis<Scalar> &); // Passed objects must be kept alive by owner
  void EmpiricalSolver();
  void addToReducedMatrix(const Eigen::Matrix<Scalar,Eigen::Dynamic,Eigen::Dynamic> &, double);
  
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

  double Mcoef_;
  
  virtual void resetSolver(int vCount, int vSize) = 0;
  virtual void assembleAndFactorReducedSystem(double Mcoef) = 0;
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
  reducedSolution_(0),
  Mcoef_(0)
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
GenDBSparsePodProjectionSolver<Scalar>::EmpiricalSolver() {
//nothing to do here
}

template <typename Scalar>
void
GenDBSparsePodProjectionSolver<Scalar>::addToReducedMatrix(const Eigen::Matrix<Scalar,Eigen::Dynamic,Eigen::Dynamic> & dummy, double dummy2){
//nothing to do
}

template <typename Scalar>
void
GenDBSparsePodProjectionSolver<Scalar>::factor() {
  for (int col = 0; col < basisSize_; ++col) {
    GenDBSparseMatrix<Scalar>::mult(projectionBasis()[col], matrixAction_[col]);
  }

  assembleAndFactorReducedSystem(Mcoef_);
}
   
template <typename Scalar>
void
GenDBSparsePodProjectionSolver<Scalar>::reSolve(GenVector<Scalar> &rhs) {
  solveReducedSystem(rhs);
}

typedef GenDBSparsePodProjectionSolver<double> DBSparsePodProjectionSolver;

} /* end namespace Rom */

#endif /* ROM_PODPROJECTIONSOLVER_H */
