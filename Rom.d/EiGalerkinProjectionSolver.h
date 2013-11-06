#ifndef ROM_EIGALERKINPROJECTIONSOLVER_H
#define ROM_EIGALERKINPROJECTIONSOLVER_H

#ifdef USE_EIGEN3
#include <Math.d/EiSparseMatrix.h>
#include "PodProjectionSolver.h"
#include "VecBasis.h"
#include <Eigen/Cholesky>
#include <iostream>

namespace Rom {

template <typename Scalar>
class GenEiSparseGalerkinProjectionSolver : public GenPodProjectionSolver<Scalar>, public GenEiSparseMatrix<Scalar> {
public:
  GenEiSparseGalerkinProjectionSolver(Connectivity *cn, DofSetArray *dsa, ConstrainedDSA *c_dsa);

  using GenEiSparseMatrix<Scalar>::neqs;

  // Full-order matrix assembly
  virtual void zeroAll();

  // Reduced-order matrix assembly
  void addReducedMass(double Mcoef);

  // Solution
  virtual void factor();
  virtual void solve(GenVector<Scalar> &rhs, GenVector<Scalar> &sol);
  virtual void reSolve(GenVector<Scalar> &rhs);

  // Reduced basis parameters
  int basisSize() const { return basisSize_; };
  const GenVecBasis<Scalar> &projectionBasis() const { return *projectionBasis_; }
  void projectionBasisIs(const GenVecBasis<Scalar> &); // Passed objects must be kept alive by owner
  
  // Data collection
  const GenVector<Scalar> &lastReducedSolution() const { 
    std::cerr << "ERROR: GenEiSparseGalerkinProjectionSolver does not implement lastReducedSolution\n";
    exit(-1);
  }
  const GenVecBasis<Scalar> &lastReducedMatrixAction() const { 
    std::cerr << "ERROR: GenEiSparseGalerkinProjectionSolver does not implement lastReducedMatrixAction\n";
    exit(-1); 
  }

private:
  int basisSize_;
  const GenVecBasis<Scalar> *projectionBasis_;
  Eigen::Matrix<Scalar,Eigen::Dynamic,Eigen::Dynamic> reducedMatrix_;
  Eigen::LLT<Eigen::Matrix<Scalar,Eigen::Dynamic,Eigen::Dynamic>, Eigen::Lower> solver_;
  
  // Disallow copy and assignment
  GenEiSparseGalerkinProjectionSolver(const GenEiSparseGalerkinProjectionSolver<Scalar> &);
  GenEiSparseGalerkinProjectionSolver<Scalar> &operator=(const GenEiSparseGalerkinProjectionSolver<Scalar> &);
};

} /* end namespace Rom */

#endif /* USE_EIGEN3 */

#endif /* ROM_EIGALERKINPROJECTIONSOLVER_H */
