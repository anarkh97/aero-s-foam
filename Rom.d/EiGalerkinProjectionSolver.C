#ifndef ROM_EIGALERKINPROJECTIONSOLVER_C
#define ROM_EIGALERKINPROJECTIONSOLVER_C

#ifdef USE_EIGEN3
#include "EiGalerkinProjectionSolver.h"
#include <stdexcept>

namespace Rom {

template <typename Scalar>
GenEiSparseGalerkinProjectionSolver<Scalar>::GenEiSparseGalerkinProjectionSolver(Connectivity *cn,
                                                   DofSetArray *dsa, ConstrainedDSA *c_dsa):
  GenEiSparseMatrix<Scalar>(cn, dsa, c_dsa), basisSize_(0), projectionBasis_(NULL)
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
GenEiSparseGalerkinProjectionSolver<Scalar>::factor()
{
  const Eigen::Matrix<Scalar,Eigen::Dynamic,Eigen::Dynamic,Eigen::ColMajor> &V = projectionBasis_->basis();
  reducedMatrix_.template triangularView<Eigen::Lower>()
    += V.transpose()*(this->M.template selfadjointView<Eigen::Upper>()*V);

  solver_.compute(reducedMatrix_);
}

template <typename Scalar>
void
GenEiSparseGalerkinProjectionSolver<Scalar>::reSolve(GenVector<Scalar> &rhs)
{
  const Eigen::Matrix<Scalar,Eigen::Dynamic,Eigen::Dynamic,Eigen::ColMajor> &V = projectionBasis_->basis();
  Eigen::Map< Eigen::Matrix<Scalar, Eigen::Dynamic, 1> > x(rhs.data(), V.rows());
  solver_.solveInPlace(x);
}

template <typename Scalar>
void
GenEiSparseGalerkinProjectionSolver<Scalar>::solve(GenVector<Scalar> &rhs, GenVector<Scalar> &sol)
{
  const Eigen::Matrix<Scalar,Eigen::Dynamic,Eigen::Dynamic,Eigen::ColMajor> &V = projectionBasis_->basis();
  Eigen::Map< Eigen::Matrix<Scalar, Eigen::Dynamic, 1> > b(rhs.data(), V.rows()), x(sol.data(), V.rows());
  x = solver_.solve(b);
}

} /* end namespace Rom */

#endif /* USE_EIGEN3 */

#endif /* ROM_EIGALERKINPROJECTIONSOLVER_C */
