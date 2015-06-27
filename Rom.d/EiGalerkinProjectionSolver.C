#ifndef ROM_EIGALERKINPROJECTIONSOLVER_C
#define ROM_EIGALERKINPROJECTIONSOLVER_C

#ifdef USE_EIGEN3
#include "EiGalerkinProjectionSolver.h"
#include <stdexcept>
#include <vector>
#include <Solvers.d/eiquadprog.hpp>
#include <Driver.d/Mpc.h>

namespace Rom {

template <typename Scalar>
GenEiSparseGalerkinProjectionSolver<Scalar>::GenEiSparseGalerkinProjectionSolver(Connectivity *cn,
                                                   DofSetArray *dsa, ConstrainedDSA *c_dsa, bool selfadjoint, double tol):
  GenEiSparseMatrix<Scalar>(cn, dsa, c_dsa, selfadjoint), 
  cdsa_(c_dsa),
  basisSize_(0), 
  dualBasisSize_(0),
  projectionBasis_(NULL), 
  dualProjectionBasis_(NULL),
  Empirical(false),
  selfadjoint_(selfadjoint),
  tol_(tol),
  startCol_(0),
  blockCols_(0)
{
}

template <typename Scalar>
void
GenEiSparseGalerkinProjectionSolver<Scalar>::setLocalBasis(int startCol, int blockCols)
{
  startCol_ = startCol;
  blockCols_ = blockCols;
  reducedMatrix_.resize(blockCols_,blockCols_);
  reducedMatrix_.setZero();
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
  reducedMatrix_ += Mcoef*Eigen::Matrix<Scalar,Eigen::Dynamic,Eigen::Dynamic>::Identity(blockCols_, blockCols_);
}

template <typename Scalar>
void
GenEiSparseGalerkinProjectionSolver<Scalar>::addToReducedMatrix(const Eigen::Matrix<Scalar,Eigen::Dynamic,Eigen::Dynamic> &ContributionMat, double Coef)
{
  reducedMatrix_ += Coef*ContributionMat;
}

template <typename Scalar>
void
GenEiSparseGalerkinProjectionSolver<Scalar>::addLMPCs(int numLMPC, LMPCons **lmpc, double Kcoef)
{
  if(numLMPC > 0 && !selfadjoint_) { std::cerr << "Error: unsymmetric solver is not supported for contact ROM\n"; exit(-1); }

  std::vector<Eigen::Triplet<Scalar> > tripletList;

  Eigen::SparseMatrix<Scalar> C(numLMPC, cdsa_->size());
  Eigen::Matrix<Scalar,Eigen::Dynamic,1> g(numLMPC);

  for(int i=0; i<numLMPC; ++i) {
    for(int j=0; j<lmpc[i]->nterms; ++j) {
      int cdof = cdsa_->locate(lmpc[i]->terms[j].nnum, 1 << lmpc[i]->terms[j].dofnum);
      if(cdof > -1) {
        tripletList.push_back(Eigen::Triplet<Scalar>(i, cdof, Scalar(lmpc[i]->terms[j].coef.r_value)));
      }
    }
    g[i] = lmpc[i]->rhs.r_value;
  }

  C.setFromTriplets(tripletList.begin(), tripletList.end());

  const Eigen::Matrix<Scalar,Eigen::Dynamic,Eigen::Dynamic,Eigen::ColMajor> &V = projectionBasis_->basis(), &W = dualProjectionBasis_->basis();
  reducedConstraintMatrix_ = Kcoef*W.transpose()*C*V;
  reducedConstraintRhs0_ = reducedConstraintRhs_ = Kcoef*W.transpose()*g;
}

template <typename Scalar>
void
GenEiSparseGalerkinProjectionSolver<Scalar>::addModalLMPCs(double Kcoef, int Wcols, std::vector<double>::const_iterator it, std::vector<double>::const_iterator it_end)
{
  std::cout << " ... Using Modal LMPCs              ..." << std::endl;
  dualBasisSize_ = Wcols;
  int counter = 0; int column = 0; int row = 0;
  reducedConstraintMatrix_.setZero(dualBasisSize_,basisSize_);
  reducedConstraintRhs0_.setZero(dualBasisSize_);
  // set reduced Constraint Matrix
  while(counter < dualBasisSize_*basisSize_) {
    reducedConstraintMatrix_(row,column) = *it; row++; counter++; it++;
    if(row == dualBasisSize_) {
      row = 0;
      column++;
    }
  }
  reducedConstraintMatrix_ *= Kcoef;
  // set reduced Constraint RHS
  row = 0;
  while(it != it_end) {
    reducedConstraintRhs0_(row) = *it; row++; it++;
  }
  reducedConstraintRhs0_ *= Kcoef;
  reducedConstraintRhs_ = reducedConstraintRhs0_;
  reducedConstraintForce_.setZero(basisSize_);
}

template <typename Scalar>
void
GenEiSparseGalerkinProjectionSolver<Scalar>::updateLMPCs(GenVector<Scalar> &_q)
{
  Eigen::Map< Eigen::Matrix<Scalar, Eigen::Dynamic, 1> > q(_q.data(), _q.size()); 
  reducedConstraintRhs_ = reducedConstraintRhs0_ - reducedConstraintMatrix_*q;
}

template <typename Scalar>
void
GenEiSparseGalerkinProjectionSolver<Scalar>::projectionBasisIs(GenVecBasis<Scalar> &reducedBasis)
{
  if (reducedBasis.vectorSize() != GenEiSparseMatrix<Scalar>::neqs()) {
    throw std::domain_error("Vectors of the reduced basis have the wrong size");
  }

  projectionBasis_ = &reducedBasis;
  basisSize_ = reducedBasis.vectorCount();
  reducedMatrix_.setZero(basisSize_, basisSize_);

  // local bases: solver uses all bases unless setLocalBasis is called
  startCol_ = 0;
  blockCols_ = basisSize_;
}

template <typename Scalar>
void
GenEiSparseGalerkinProjectionSolver<Scalar>::dualProjectionBasisIs(GenVecBasis<Scalar> &dualReducedBasis)
{
  dualProjectionBasis_ = &dualReducedBasis;
  dualBasisSize_ = dualReducedBasis.vectorCount();
  reducedConstraintMatrix_.setZero(dualBasisSize_, basisSize_);
  reducedConstraintForce_.setZero(basisSize_);
}

template <typename Scalar>
void
GenEiSparseGalerkinProjectionSolver<Scalar>::EmpiricalSolver()
{
  std::cout << " ... Empirical Solver Selected      ..." << std::endl;
  Empirical = true;
}

template <typename Scalar>
void
GenEiSparseGalerkinProjectionSolver<Scalar>::factor()
{
  LocalBasisType V = projectionBasis_->basis().block(0,startCol_,projectionBasis_->size(),blockCols_);

  if(selfadjoint_ && !Empirical) {
    reducedMatrix_.template triangularView<Eigen::Lower>()
    += V.transpose()*(this->M.template selfadjointView<Eigen::Upper>()*V);
    if(dualBasisSize_ > 0) c1_ = reducedMatrix_.trace();
    llt_.compute(reducedMatrix_);
  }
  else {
    if(Empirical) {
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
  LocalBasisType V = projectionBasis_->basis().block(0,startCol_,projectionBasis_->size(),blockCols_);
  Eigen::Map< Eigen::Matrix<Scalar, Eigen::Dynamic, 1> > x(rhs.data()+startCol_, V.cols());

  if(dualBasisSize_ > 0) {
    Eigen::Matrix<Scalar,Eigen::Dynamic,Eigen::Dynamic> CE(0,0), CI = -reducedConstraintMatrix_.block(0,startCol_,dualBasisSize_,blockCols_).transpose();
    Eigen::Matrix<Scalar,Eigen::Dynamic,1> g0 = -x, ce0(0,1), _x(V.cols()), Lambda(0,1), Mu(dualBasisSize_);
    solve_quadprog2(llt_, c1_, g0, CE, ce0, CI, reducedConstraintRhs_, _x, &Lambda, &Mu, tol_);
    x = _x;
    reducedConstraintForce_.setZero();
    reducedConstraintForce_.segment(startCol_,blockCols_) = reducedConstraintMatrix_.block(0,startCol_,dualBasisSize_,blockCols_).transpose()*Mu;
  }
  else if(selfadjoint_ && !Empirical) llt_.solveInPlace(x);
  else x = (lu_.solve(x)).eval();

  for(int i=0; i<startCol_; ++i) rhs[i] = 0;
  for(int i=startCol_+blockCols_; i<rhs.size(); ++i) rhs[i] = 0;
}

template <typename Scalar>
void
GenEiSparseGalerkinProjectionSolver<Scalar>::solve(GenVector<Scalar> &rhs, GenVector<Scalar> &sol)
{
  LocalBasisType V = projectionBasis_->basis().block(0,startCol_,projectionBasis_->size(),blockCols_);
  Eigen::Map< Eigen::Matrix<Scalar, Eigen::Dynamic, 1> > b(rhs.data()+startCol_, V.cols()), x(sol.data()+startCol_, V.cols());
  sol.zero();

  if(dualBasisSize_ > 0) {
    Eigen::Matrix<Scalar,Eigen::Dynamic,Eigen::Dynamic> CE(0,0), CI = -reducedConstraintMatrix_.block(0,startCol_,dualBasisSize_,blockCols_).transpose();
    Eigen::Matrix<Scalar,Eigen::Dynamic,1> g0 = -b, ce0(0,1), _x(V.cols()), Lambda(0,1), Mu(dualBasisSize_);
    solve_quadprog2(llt_, c1_, g0, CE, ce0, CI, reducedConstraintRhs_, _x, &Lambda, &Mu, tol_);
    x = _x;
    reducedConstraintForce_.setZero();
    reducedConstraintForce_.segment(startCol_,blockCols_) = reducedConstraintMatrix_.block(0,startCol_,dualBasisSize_,blockCols_).transpose()*Mu;
  }
  else if(selfadjoint_ && !Empirical) x = llt_.solve(b);
  else x = lu_.solve(b);
}

template <typename Scalar>
double
GenEiSparseGalerkinProjectionSolver<Scalar>::getResidualNorm(const GenVector<Scalar> &v)
{
  Eigen::Map< Eigen::Matrix<Scalar, Eigen::Dynamic, 1> > vj(v.data()+startCol_, blockCols_);
  return vj.norm();
}

} /* end namespace Rom */

#endif /* USE_EIGEN3 */

#endif /* ROM_EIGALERKINPROJECTIONSOLVER_C */
