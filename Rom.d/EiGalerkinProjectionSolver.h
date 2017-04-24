#ifndef ROM_EIGALERKINPROJECTIONSOLVER_H
#define ROM_EIGALERKINPROJECTIONSOLVER_H

#ifdef USE_EIGEN3
#include <Math.d/EiSparseMatrix.h>
#include <Timers.d/Timing.h>
#include "PodProjectionSolver.h"
#include "VecBasis.h"
#include "DistrVecBasis.h"
#include <Eigen/Dense>
#include <iostream>

template <class Scalar> class GenSubDOp;

namespace Rom {

template <typename Scalar, template<typename> class GenVecType = GenVector, class BaseSolver = GenPodProjectionSolver<Scalar, GenVecType> >
class GenEiSparseGalerkinProjectionSolver : public BaseSolver, public GenEiSparseMatrix<Scalar>
{
  typedef const Eigen::Block<const Eigen::Matrix<Scalar,Eigen::Dynamic,Eigen::Dynamic,Eigen::ColMajor>, 
                             Eigen::Dynamic, Eigen::Dynamic, false> LocalBasisType;
public:
  GenEiSparseGalerkinProjectionSolver(Connectivity *cn, DofSetArray *dsa, ConstrainedDSA *c_dsa, bool = true, double tol = 1e-6);
  GenEiSparseGalerkinProjectionSolver(Connectivity *cn, DofSetArray *dsa, ConstrainedDSA *c_dsa, int numSub_, GenSparseMatrix<Scalar>**, bool = true, double tol = 1e-6);
  ~GenEiSparseGalerkinProjectionSolver() { /* do nothing */ }

  using GenEiSparseMatrix<Scalar>::neqs;

  // Local bases
  void setLocalBasis(int startCol, int blockCols);
  void setLocalDualBasis(int startCol, int blockCols);

  // Full-order matrix assembly
  virtual void zeroAll();

  // member function for GenParallelSolver inheritance (sp?)
  double getSolutionTime() { return solveTime; }
  Timings& getTimers() { return times; }

  // Reduced-order matrix assembly
  void addReducedMass(double Mcoef);

  // Constraint assembly
  void activateContact() { contact_ = true; }
  void addLMPCs(int numLMPC, LMPCons **lmpc, double Kcoef);
  void addModalLMPCs(double Kcoef, int Wcols,std::vector<double>::const_iterator it, std::vector<double>::const_iterator it_end);
  void updateLMPCs(GenVecType<Scalar> &q);

  // Solution
  virtual void factor();
  virtual void refactor();
  virtual void solve(GenVecType<Scalar> &rhs, GenVecType<Scalar> &sol);
  virtual void reSolve(GenVecType<Scalar> &rhs);

  /*double refactorTime = 0.0;
  double reducedSolveTime = 0.0;
  int nRefactor = 0; 
  int nSolve = 0;*/

  // Reduced basis parameters
  int basisSize() const { return basisSize_; };
  GenVecBasis<Scalar,GenVecType> &projectionBasis() { return *projectionBasis_; }
  GenVecBasis<Scalar,GenVecType> &dualProjectionBasis() { return *dualProjectionBasis_; }
  void projectionBasisIs(GenVecBasis<Scalar,GenVecType> &); // Passed objects must be kept alive by owner
  void dualProjectionBasisIs(GenVecBasis<Scalar,GenVecType> &);
  void storeReducedMass(Eigen::Matrix<Scalar,Eigen::Dynamic,Eigen::Dynamic> &);
  void EmpiricalSolver(); 
  void addToReducedMatrix(const Eigen::Matrix<Scalar,Eigen::Dynamic,Eigen::Dynamic> &, double = 1.0); 

  // Data collection
  const GenVecType<Scalar> &lastReducedSolution() const { 
    std::cerr << "ERROR: GenEiSparseGalerkinProjectionSolver does not implement lastReducedSolution\n";
    exit(-1);
  }
  const GenVecBasis<Scalar,GenVecType> &lastReducedMatrixAction() const { 
    std::cerr << "ERROR: GenEiSparseGalerkinProjectionSolver does not implement lastReducedMatrixAction\n";
    exit(-1); 
  }
  const Eigen::Matrix<Scalar,Eigen::Dynamic,1> &lastReducedConstraintForce() const {
    return reducedConstraintForce_;
  }
  double getResidualNorm(const GenVecType<Scalar> &v);
  double getFNormSq(GenVecType<Scalar> &v);

  GenSparseMatrix<Scalar> * getSpMat(int i) { return spMat[i]; }

private:
  ConstrainedDSA *cdsa_;
  GenSparseMatrix<Scalar> **spMat;
  GenSubDOp<Scalar> *K;
  int numSub; 
  bool selfadjoint_;
  bool Empirical;
  bool contact_;
  int basisSize_, dualBasisSize_; // global basis quantities
  int startCol_, blockCols_; // local bases quantities
  int startDualCol_, dualBlockCols_;
  double tol_; // convergence tolerance used by QP solver for contact
  //double solveTime;
  //Timings times; 
  Scalar c1_; // trace of reducedConstraintMatrix_
  GenVecBasis<Scalar,GenVecType> *projectionBasis_, *dualProjectionBasis_;
  Eigen::Matrix<Scalar,Eigen::Dynamic,Eigen::Dynamic> reducedMatrix_, reducedConstraintMatrix_;
  Eigen::Matrix<Scalar,Eigen::Dynamic,1> reducedConstraintRhs_, reducedConstraintRhs0_, reducedConstraintForce_, VtMV;
  Eigen::LLT<Eigen::Matrix<Scalar,Eigen::Dynamic,Eigen::Dynamic>, Eigen::Lower> llt_;
  Eigen::PartialPivLU<Eigen::Matrix<Scalar,Eigen::Dynamic,Eigen::Dynamic> > lu_;
  
  // Disallow copy and assignment
  GenEiSparseGalerkinProjectionSolver(const GenEiSparseGalerkinProjectionSolver<Scalar> &);
  GenEiSparseGalerkinProjectionSolver<Scalar,GenVecType> &operator=(const GenEiSparseGalerkinProjectionSolver<Scalar> &);
};

} /* end namespace Rom */

#endif /* USE_EIGEN3 */

#endif /* ROM_EIGALERKINPROJECTIONSOLVER_H */
