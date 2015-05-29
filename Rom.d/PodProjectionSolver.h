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

  // Constraint assembly
  virtual void addLMPCs(int numLMPC, LMPCons **lmpc, double Kcoef) = 0;
  virtual void addModalLMPCs(double Kcoef, int Wcols, std::vector<double>::const_iterator it, std::vector<double>::const_iterator it_end) = 0;
  virtual void updateLMPCs(GenVector<Scalar> &q) = 0;

  // Solution
  virtual void factor() = 0;
  virtual void reSolve(GenVector<Scalar> &rhs) = 0;

  // Reduced basis parameters
  virtual int basisSize() const = 0;
  virtual const GenVecBasis<Scalar> &projectionBasis() const = 0;
  virtual void projectionBasisIs(const GenVecBasis<Scalar> &) = 0; 
  virtual void dualProjectionBasisIs(const GenVecBasis<Scalar> &) = 0;
  virtual void EmpiricalSolver() = 0; 
#ifdef USE_EIGEN3
  virtual void addToReducedMatrix(const Eigen::Matrix<Scalar,Eigen::Dynamic,Eigen::Dynamic> &, double = 1.0) = 0;
#endif

  // Data collection
  virtual const GenVector<Scalar> &lastReducedSolution() const = 0;
  virtual const GenVecBasis<Scalar> &lastReducedMatrixAction() const = 0;
#ifdef USE_EIGEN3
  virtual const Eigen::Matrix<Scalar,Eigen::Dynamic,1> &lastReducedConstraintForce() const = 0;
#endif
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

  // Constraint assembly
  void addLMPCs(int numLMPC, LMPCons **lmpc, double Kcoef) {}
  void addModalLMPCs(double Kcoef, int Wcols,std::vector<double>::const_iterator it, std::vector<double>::const_iterator it_end) {}
  void updateLMPCs(GenVector<Scalar> &q) {}

  // Solution
  virtual void factor();
  virtual void reSolve(GenVector<Scalar> &rhs);

  // Reduced basis parameters
  int basisSize() const { return basisSize_; }
  const GenVecBasis<Scalar> &projectionBasis() const { return *projectionBasis_; }
  void projectionBasisIs(const GenVecBasis<Scalar> &); // Passed objects must be kept alive by owner
  void dualProjectionBasisIs(const GenVecBasis<Scalar> &);
  void EmpiricalSolver();
#ifdef USE_EIGEN3
  void addToReducedMatrix(const Eigen::Matrix<Scalar,Eigen::Dynamic,Eigen::Dynamic> &, double);
#endif
  
  // Data collection
  const GenVector<Scalar> &lastReducedSolution() const { return reducedSolution_; }
  const GenVecBasis<Scalar> &lastReducedMatrixAction() const { return matrixAction_; }
#ifdef USE_EIGEN3
  const Eigen::Matrix<Scalar,Eigen::Dynamic,1> &lastReducedConstraintForce() const {
    std::cerr << "ERROR: GenDBSparsePodProjectionSolver does not implement lastReducedConstraintForce\n";
    exit(-1);
  }
#endif

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
GenDBSparsePodProjectionSolver<Scalar>::dualProjectionBasisIs(const GenVecBasis<Scalar> &) {
  throw std::domain_error("Selected solver does not use dual projection basis\n");
}

template <typename Scalar>
void
GenDBSparsePodProjectionSolver<Scalar>::EmpiricalSolver() {
//nothing to do here
}

#ifdef USE_EIGEN3
template <typename Scalar>
void
GenDBSparsePodProjectionSolver<Scalar>::addToReducedMatrix(const Eigen::Matrix<Scalar,Eigen::Dynamic,Eigen::Dynamic> & dummy, double dummy2){
//nothing to do
}
#endif

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
