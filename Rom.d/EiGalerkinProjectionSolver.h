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
  typedef const Eigen::Block<const Eigen::Matrix<Scalar,Eigen::Dynamic,Eigen::Dynamic,Eigen::ColMajor>, 
                             Eigen::Dynamic, Eigen::Dynamic, false> LocalBasisType;
public:
  GenEiSparseGalerkinProjectionSolver(Connectivity *cn, DofSetArray *dsa, ConstrainedDSA *c_dsa, bool = true, double tol = 1e-6);

  using GenEiSparseMatrix<Scalar>::neqs;

  // Local bases
  void setLocalBasis(int startCol, int blockCols);
  void setLocalDualBasis(int startCol, int blockCols);

  // Full-order matrix assembly
  virtual void zeroAll();

  // Reduced-order matrix assembly
  void addReducedMass(double Mcoef);

  // Constraint assembly
  void activateContact() { contact_ = true; }
  void addLMPCs(int numLMPC, LMPCons **lmpc, double Kcoef);
  void addModalLMPCs(double Kcoef, int Wcols,std::vector<double>::const_iterator it, std::vector<double>::const_iterator it_end);
  void updateLMPCs(GenVector<Scalar> &q);

  // Solution
  virtual void factor();
  virtual void solve(GenVector<Scalar> &rhs, GenVector<Scalar> &sol);
  virtual void reSolve(GenVector<Scalar> &rhs);

  // Reduced basis parameters
  int basisSize() const { return basisSize_; };
  GenVecBasis<Scalar> &projectionBasis() { return *projectionBasis_; }
  GenVecBasis<Scalar> &dualProjectionBasis() { return *dualProjectionBasis_; }
  void projectionBasisIs(GenVecBasis<Scalar> &); // Passed objects must be kept alive by owner
  void dualProjectionBasisIs(GenVecBasis<Scalar> &);
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
  double getResidualNorm(const GenVector<Scalar> &v);

private:
  ConstrainedDSA *cdsa_;
  bool selfadjoint_;
  bool Empirical;
  bool contact_;
  int basisSize_, dualBasisSize_; // global basis quantities
  int startCol_, blockCols_; // local bases quantities
  int startDualCol_, dualBlockCols_;
  double tol_; // convergence tolerance used by QP solver for contact
  Scalar c1_; // trace of reducedConstraintMatrix_
  GenVecBasis<Scalar> *projectionBasis_, *dualProjectionBasis_;
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
