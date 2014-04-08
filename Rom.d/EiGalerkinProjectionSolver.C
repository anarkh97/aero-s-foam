#ifndef ROM_EIGALERKINPROJECTIONSOLVER_C
#define ROM_EIGALERKINPROJECTIONSOLVER_C

#ifdef USE_EIGEN3
#include "EiGalerkinProjectionSolver.h"
#include <stdexcept>

namespace Rom {

template <typename Scalar>
GenEiSparseGalerkinProjectionSolver<Scalar>::GenEiSparseGalerkinProjectionSolver(Connectivity *cn,
                                                   DofSetArray *dsa, ConstrainedDSA *c_dsa, bool selfadjoint):
  GenEiSparseMatrix<Scalar>(cn, dsa, c_dsa, selfadjoint), 
  basisSize_(0), 
  projectionBasis_(NULL), 
  Empirical(false),
  selfadjoint_(selfadjoint)
{
}

template <typename Scalar>
void
GenEiSparseGalerkinProjectionSolver<Scalar>::zeroAll()
{
  GenEiSparseMatrix<Scalar>::zeroAll();
  reducedMatrix_.setZero();
}

template <typename Scalar>
void
GenEiSparseGalerkinProjectionSolver<Scalar>::addReducedMass(double Mcoef)
{
  reducedMatrix_ += Mcoef*Eigen::Matrix<Scalar,Eigen::Dynamic,Eigen::Dynamic>::Identity(basisSize_, basisSize_);
}

template <typename Scalar>
void
GenEiSparseGalerkinProjectionSolver<Scalar>::addToReducedMatrix(const Eigen::Matrix<Scalar,Eigen::Dynamic,Eigen::Dynamic> &ContributionMat, double Coef){

  reducedMatrix_ += Coef*ContributionMat;

}

template <typename Scalar>
void
GenEiSparseGalerkinProjectionSolver<Scalar>::projectionBasisIs(const GenVecBasis<Scalar> &reducedBasis)
{
  if (reducedBasis.vectorSize() != GenEiSparseMatrix<Scalar>::neqs()) {
    throw std::domain_error("Vectors of the reduced basis have the wrong size");
  }

  projectionBasis_ = &reducedBasis;
  basisSize_ = reducedBasis.vectorCount();
  reducedMatrix_.setZero(basisSize_, basisSize_);
}

template <typename Scalar>
void
GenEiSparseGalerkinProjectionSolver<Scalar>::EmpiricalSolver()
{
  std::cout << " ... Empirical Solver Selected ..." << std::endl;
  Empirical = true;
}

template <typename Scalar>
void
GenEiSparseGalerkinProjectionSolver<Scalar>::factor()
{
  const Eigen::Matrix<Scalar,Eigen::Dynamic,Eigen::Dynamic,Eigen::ColMajor> &V = projectionBasis_->basis();
  if(selfadjoint_ && !Empirical) {
    if(Empirical){
      llt_.compute(reducedMatrix_);//this is unstable for Empirical methods, symmetry is lost
    }else {
      reducedMatrix_.template triangularView<Eigen::Lower>()
      += V.transpose()*(this->M.template selfadjointView<Eigen::Upper>()*V);
      llt_.compute(reducedMatrix_);
    }
  } else {
    if(Empirical){
      lu_.compute(reducedMatrix_);
    } else {
      reducedMatrix_ += V.transpose()*(this->M*V);
      lu_.compute(reducedMatrix_);
    }
  }
}

template <typename Scalar>
void
GenEiSparseGalerkinProjectionSolver<Scalar>::reSolve(GenVector<Scalar> &rhs)
{
  const Eigen::Matrix<Scalar,Eigen::Dynamic,Eigen::Dynamic,Eigen::ColMajor> &V = projectionBasis_->basis();
  Eigen::Map< Eigen::Matrix<Scalar, Eigen::Dynamic, 1> > x(rhs.data(), V.cols());
  if(selfadjoint_ && !Empirical) llt_.solveInPlace(x);
  else x = (lu_.solve(x)).eval();
}

template <typename Scalar>
void
GenEiSparseGalerkinProjectionSolver<Scalar>::solve(GenVector<Scalar> &rhs, GenVector<Scalar> &sol)
{
  const Eigen::Matrix<Scalar,Eigen::Dynamic,Eigen::Dynamic,Eigen::ColMajor> &V = projectionBasis_->basis();
  Eigen::Map< Eigen::Matrix<Scalar, Eigen::Dynamic, 1> > b(rhs.data(), V.cols()), x(sol.data(), V.cols());
  if(selfadjoint_ && !Empirical) x = llt_.solve(b);
  else x = lu_.solve(b);
}

} /* end namespace Rom */

#endif /* USE_EIGEN3 */

#endif /* ROM_EIGALERKINPROJECTIONSOLVER_C */
