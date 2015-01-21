#ifndef ROM_EIGALERKINPROJECTIONSOLVER_H
#define ROM_EIGALERKINPROJECTIONSOLVER_H

#ifdef USE_EIGEN3
#include <Math.d/EiSparseMatrix.h>
#include "PodProjectionSolver.h"
#include "VecBasis.h"
#include <Eigen/Dense>
#include <iostream>

namespace Rom {

template <typename Scalar>
class GenEiSparseGalerkinProjectionSolver : public GenPodProjectionSolver<Scalar>, public GenEiSparseMatrix<Scalar> {
public:
  GenEiSparseGalerkinProjectionSolver(Connectivity *cn, DofSetArray *dsa, ConstrainedDSA *c_dsa, bool = true, double tol = 1e-6);

  using GenEiSparseMatrix<Scalar>::neqs;

  // Full-order matrix assembly
  virtual void zeroAll();

  // Reduced-order matrix assembly
  void addReducedMass(double Mcoef);

  // Constraint assembly
  void addLMPCs(int numLMPC, LMPCons **lmpc, double Kcoef);
  void updateLMPCs(GenVector<Scalar> &q);

  // Solution
  virtual void factor();
  virtual void solve(GenVector<Scalar> &rhs, GenVector<Scalar> &sol);
  virtual void reSolve(GenVector<Scalar> &rhs);

  // Reduced basis parameters
  int basisSize() const { return basisSize_; };
  const GenVecBasis<Scalar> &projectionBasis() const { return *projectionBasis_; }
  void projectionBasisIs(const GenVecBasis<Scalar> &); // Passed objects must be kept alive by owner
  void dualProjectionBasisIs(const GenVecBasis<Scalar> &);
  void EmpiricalSolver(); 
  void addToReducedMatrix(const Eigen::Matrix<Scalar,Eigen::Dynamic,Eigen::Dynamic> &, double = 1.0); 

  // Data collection
  const GenVector<Scalar> &lastReducedSolution() const { 
    std::cerr << "ERROR: GenEiSparseGalerkinProjectionSolver does not implement lastReducedSolution\n";
    exit(-1);
  }
  const GenVecBasis<Scalar> &lastReducedMatrixAction() const { 
    std::cerr << "ERROR: GenEiSparseGalerkinProjectionSolver does not implement lastReducedMatrixAction\n";
    exit(-1); 
  }
  const Eigen::Matrix<Scalar,Eigen::Dynamic,1> &lastReducedConstraintForce() const {
    return reducedConstraintForce_;
  }

private:
  ConstrainedDSA *cdsa_;
  bool selfadjoint_;
  bool Empirical;
  int basisSize_, dualBasisSize_;
  double tol_; // convergence tolerance used by QP solver for contact
  Scalar c1_; // trace of reducedConstraintMatrix_
  const GenVecBasis<Scalar> *projectionBasis_, *dualProjectionBasis_;
  Eigen::Matrix<Scalar,Eigen::Dynamic,Eigen::Dynamic> reducedMatrix_, reducedConstraintMatrix_;
  Eigen::Matrix<Scalar,Eigen::Dynamic,1> reducedConstraintRhs_, reducedConstraintRhs0_, reducedConstraintForce_;
  Eigen::LLT<Eigen::Matrix<Scalar,Eigen::Dynamic,Eigen::Dynamic>, Eigen::Lower> llt_;
  Eigen::PartialPivLU<Eigen::Matrix<Scalar,Eigen::Dynamic,Eigen::Dynamic> > lu_;
  
  // Disallow copy and assignment
  GenEiSparseGalerkinProjectionSolver(const GenEiSparseGalerkinProjectionSolver<Scalar> &);
  GenEiSparseGalerkinProjectionSolver<Scalar> &operator=(const GenEiSparseGalerkinProjectionSolver<Scalar> &);
};

} /* end namespace Rom */

#endif /* USE_EIGEN3 */

#endif /* ROM_EIGALERKINPROJECTIONSOLVER_H */
