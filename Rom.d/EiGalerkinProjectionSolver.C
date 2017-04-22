#ifndef ROM_EIGALERKINPROJECTIONSOLVER_C
#define ROM_EIGALERKINPROJECTIONSOLVER_C

#ifdef USE_EIGEN3
#include "EiGalerkinProjectionSolver.h"
#include "DistrVecBasisOps.h"
#include "PodProjectionSolver.h"
#include <stdexcept>
#include <vector>
#include <Solvers.d/eiquadprog.hpp>
#include <Solvers.d/ParallelSolver.h>
#include <Driver.d/Mpc.h>
#include <Paral.d/SubDOp.h>

namespace Rom {

template <typename Scalar, template<typename> class GenVecType, class BaseSolver>
GenEiSparseGalerkinProjectionSolver<Scalar,GenVecType,BaseSolver>::GenEiSparseGalerkinProjectionSolver(Connectivity *cn,
                                                   DofSetArray *dsa, ConstrainedDSA *c_dsa, bool selfadjoint, double tol):
  GenEiSparseMatrix<Scalar>(cn, dsa, c_dsa, selfadjoint), 
  spMat(NULL),
  K(NULL),
  cdsa_(c_dsa),
  basisSize_(0), 
  dualBasisSize_(0),
  projectionBasis_(NULL), 
  dualProjectionBasis_(NULL),
  Empirical(false),
  selfadjoint_(selfadjoint),
  contact_(false),
  tol_(tol),
  solveTime(0.0),
  startCol_(0),
  blockCols_(0),
  startDualCol_(0),
  dualBlockCols_(0)
{
}

template <typename Scalar, template<typename> class GenVecType, class BaseSolver>
GenEiSparseGalerkinProjectionSolver<Scalar,GenVecType,BaseSolver>::GenEiSparseGalerkinProjectionSolver(Connectivity *cn,
                                                           DofSetArray *dsa, ConstrainedDSA *c_dsa, int numSub_,
                                                           GenSparseMatrix<Scalar> **spMat_, bool selfadjoint, double tol):
  GenEiSparseMatrix<Scalar>(cn, dsa, c_dsa, selfadjoint),
  spMat(spMat_),
  cdsa_(c_dsa),
  basisSize_(0),
  numSub(numSub_),
  dualBasisSize_(0),
  projectionBasis_(NULL), 
  dualProjectionBasis_(NULL),
  Empirical(false),
  selfadjoint_(selfadjoint),
  contact_(false),
  tol_(tol),
  solveTime(0.0),
  startCol_(0),
  blockCols_(0),
  startDualCol_(0),
  dualBlockCols_(0)
{
  K = new GenSubDOp<Scalar>(numSub,spMat);
}

template <typename Scalar, template<typename> class GenVecType, class BaseSolver>
void
GenEiSparseGalerkinProjectionSolver<Scalar,GenVecType,BaseSolver>::setLocalBasis(int startCol, int blockCols) 
{
  startCol_  = startCol;
  blockCols_ = blockCols;
  reducedMatrix_.resize(blockCols_,blockCols_);
  reducedMatrix_.setZero();
}

template <typename Scalar, template<typename> class GenVecType, class BaseSolver>
void
GenEiSparseGalerkinProjectionSolver<Scalar,GenVecType,BaseSolver>::setLocalDualBasis(int startDualCol, int dualBlockCols)
{
  startDualCol_ = startDualCol;  // set which columns are to be used in the reduced constraint matrix
  dualBlockCols_ = dualBlockCols;
}

template <typename Scalar, template<typename> class GenVecType, class BaseSolver>
void
GenEiSparseGalerkinProjectionSolver<Scalar,GenVecType,BaseSolver>::zeroAll()
{
  GenEiSparseMatrix<Scalar>::zeroAll();
  reducedMatrix_.setZero();
  if(K) K->zeroAll();
}

template <typename Scalar, template<typename> class GenVecType, class BaseSolver>
void
GenEiSparseGalerkinProjectionSolver<Scalar,GenVecType,BaseSolver>::storeReducedMass(Eigen::Matrix<Scalar,Eigen::Dynamic,Eigen::Dynamic> &VtMV_)
{
  VtMV = VtMV_;
}

template <typename Scalar, template<typename> class GenVecType, class BaseSolver>
void
GenEiSparseGalerkinProjectionSolver<Scalar,GenVecType,BaseSolver>::addReducedMass(double Mcoef)
{
  if(VtMV.rows() > 0)
    reducedMatrix_ += Mcoef*VtMV.block(startCol_,startCol_,blockCols_,blockCols_);
  else
    reducedMatrix_ += Mcoef*Eigen::Matrix<Scalar,Eigen::Dynamic,Eigen::Dynamic>::Identity(blockCols_, blockCols_);
}

template <typename Scalar, template<typename> class GenVecType, class BaseSolver>
void
GenEiSparseGalerkinProjectionSolver<Scalar,GenVecType,BaseSolver>::addToReducedMatrix(const Eigen::Matrix<Scalar,Eigen::Dynamic,Eigen::Dynamic> &ContributionMat, double Coef)
{
  reducedMatrix_ += Coef*ContributionMat.block(startCol_,startCol_,blockCols_,blockCols_);
}

template <typename Scalar, template<typename> class GenVecType, class BaseSolver>
void
GenEiSparseGalerkinProjectionSolver<Scalar,GenVecType,BaseSolver>::addLMPCs(int numLMPC, LMPCons **lmpc, double Kcoef)
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

  // this is only called once, so set W and V to all columns, if using local basis for either W or V, then select the correct block diagonal element
  const Eigen::Matrix<Scalar,Eigen::Dynamic,Eigen::Dynamic,Eigen::ColMajor> &V = projectionBasis_->basis();
  const Eigen::Matrix<Scalar,Eigen::Dynamic,Eigen::Dynamic,Eigen::ColMajor> &W = dualProjectionBasis_->basis();
  reducedConstraintMatrix_ = Kcoef*W.transpose()*C*V; 
  reducedConstraintRhs0_ = reducedConstraintRhs_ = Kcoef*W.transpose()*g; 
}

template <typename Scalar, template<typename> class GenVecType, class BaseSolver>
void
GenEiSparseGalerkinProjectionSolver<Scalar,GenVecType,BaseSolver>::addModalLMPCs(double Kcoef, int Wcols, std::vector<double>::const_iterator it, std::vector<double>::const_iterator it_end)
{
  //std::cout << " ... Using Modal LMPCs              ..." << std::endl;
  filePrint(stderr," ... Using Modal LMPCs              ...\n");
  dualBasisSize_  = dualBlockCols_ = Wcols;
  int counter = 0; int column = 0; int row = 0;
  reducedConstraintMatrix_.setZero(dualBasisSize_,basisSize_); //allocate enough space for all local bases
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

template <typename Scalar, template<typename> class GenVecType, class BaseSolver>
void
GenEiSparseGalerkinProjectionSolver<Scalar,GenVecType,BaseSolver>::updateLMPCs(GenVecType<Scalar> &_q)
{
  Eigen::Map< Eigen::Matrix<Scalar, Eigen::Dynamic, 1> > q(_q.data(), _q.size()); 
  reducedConstraintRhs_ = reducedConstraintRhs0_ - reducedConstraintMatrix_*q;
}

template <typename Scalar, template<typename> class GenVecType, class BaseSolver>
void
GenEiSparseGalerkinProjectionSolver<Scalar,GenVecType,BaseSolver>::projectionBasisIs(GenVecBasis<Scalar,GenVecType> &reducedBasis)
{
/*  if (reducedBasis.globalVectorSize() != GenEiSparseMatrix<Scalar>::neqs()) {
    fprintf(stderr,"Basis Vector length: %d, Number of Unconstrained DOFs: %d \n",reducedBasis.globalVectorSize(), GenEiSparseMatrix<Scalar>::neqs());
    throw std::domain_error("Vectors of the reduced basis have the wrong size");
  }*/

  projectionBasis_ = &reducedBasis;
  basisSize_ = reducedBasis.vectorCount();
  reducedMatrix_.setZero(basisSize_, basisSize_);

  // local bases: solver uses all bases unless setLocalBasis is called
  startCol_ = 0;
  blockCols_ = basisSize_;
}

template <typename Scalar, template<typename> class GenVecType, class BaseSolver>
void
GenEiSparseGalerkinProjectionSolver<Scalar,GenVecType,BaseSolver>::dualProjectionBasisIs(GenVecBasis<Scalar,GenVecType> &dualReducedBasis)
{
  dualProjectionBasis_ = &dualReducedBasis;
  dualBasisSize_       = dualReducedBasis.vectorCount();
  reducedConstraintMatrix_.setZero(dualBasisSize_, basisSize_);
  reducedConstraintForce_.setZero(basisSize_);

  // local basis: dual solver uses all columns unless setLocalDualBasis is called
  startDualCol_  = 0; 
  dualBlockCols_ = dualBasisSize_;
}

template <typename Scalar, template<typename> class GenVecType, class BaseSolver>
void
GenEiSparseGalerkinProjectionSolver<Scalar,GenVecType,BaseSolver>::EmpiricalSolver()
{
  std::cout << " ... Empirical Solver Selected      ..." << std::endl;
  Empirical = true;
}

template <typename Scalar, template<typename> class GenVecType, class BaseSolver>
void
GenEiSparseGalerkinProjectionSolver<Scalar,GenVecType,BaseSolver>::factor()
{
  LocalBasisType V = projectionBasis_->basis().block(0,startCol_,projectionBasis_->size(),blockCols_);

  if(selfadjoint_ && !Empirical) {
    reducedMatrix_.template triangularView<Eigen::Lower>()
    += V.transpose()*(this->M.template selfadjointView<Eigen::Upper>()*V);
    if(dualBlockCols_ > 0) c1_ = reducedMatrix_.trace();
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

template <typename Scalar, template<typename> class GenVecType, class BaseSolver>
void
GenEiSparseGalerkinProjectionSolver<Scalar,GenVecType,BaseSolver>::refactor()
{
}

template <>
void
GenEiSparseGalerkinProjectionSolver<double,GenDistrVector,GenParallelSolver<double> >::refactor()
{ // for parallel implicit ROM, V^T*M*V is computed in problem descriptor and passed through
  // addReducedMass and addToReducedMatrix
  double dummyTime = 0.0;
  if(selfadjoint_ && !Empirical) {
   
    GenFullSquareMatrix<double> K_reduced; // local data structure
    calculateReducedStiffness(*K, *projectionBasis_, K_reduced, selfadjoint_); // parallel multiplication

    // upper half of K_reduced is filled, but data structure is row major, Krmap is Column major
    Eigen::Map< Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> > Krmap(K_reduced.data(), blockCols_, blockCols_);
    reducedMatrix_.triangularView<Eigen::Lower>() 
    += Krmap; 

    dummyTime -= getTime();
    if(dualBlockCols_ > 0) c1_ = reducedMatrix_.trace();
    llt_.compute(reducedMatrix_);
    dummyTime += getTime();
  }
  else {
    if(Empirical) {
      lu_.compute(reducedMatrix_); 
    } else { 

      GenFullSquareMatrix<double> K_reduced;
      calculateReducedStiffness(*K, *projectionBasis_, K_reduced);

      Eigen::Map< Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> > Krmap(K_reduced.data(), blockCols_, blockCols_);
      reducedMatrix_ += Krmap;

      dummyTime -= getTime(); 
      lu_.compute(reducedMatrix_);
      dummyTime += getTime();
    }
  }
  refactorTime += dummyTime/1000.0;
  nRefactor++;
  filePrint(stderr,"Average Refactor Time: %3.2e \n",refactorTime/double(nRefactor)); 
}

template <typename Scalar, template<typename> class GenVecType, class BaseSolver>
void
GenEiSparseGalerkinProjectionSolver<Scalar,GenVecType,BaseSolver>::reSolve(GenVecType<Scalar> &rhs)
{
  Eigen::Map< Eigen::Matrix<Scalar, Eigen::Dynamic, 1> > x(rhs.data()+startCol_, blockCols_);

  double dummyTime = -1.0*getTime();
  if(dualBlockCols_ > 0 && contact_) {
    Eigen::Matrix<Scalar,Eigen::Dynamic,Eigen::Dynamic> CE(0,0), CI = -reducedConstraintMatrix_.block(startDualCol_,startCol_,dualBlockCols_,blockCols_).transpose();
    Eigen::Matrix<Scalar,Eigen::Dynamic,1> g0 = -x, ce0(0,1), _x(blockCols_), Lambda(0,1), Mu(dualBlockCols_);
    solve_quadprog2(llt_, c1_, g0, CE, ce0, CI, reducedConstraintRhs_.segment(startDualCol_,dualBlockCols_), _x, &Lambda, &Mu, tol_);
    x = _x;
    reducedConstraintForce_.setZero();
    reducedConstraintForce_.segment(startCol_,blockCols_) = reducedConstraintMatrix_.block(startDualCol_,startCol_,dualBlockCols_,blockCols_).transpose()*Mu;
  }
  else if(selfadjoint_ && !Empirical) llt_.solveInPlace(x);
  else x = (lu_.solve(x)).eval();
  
  dummyTime += getTime();
  reducedSolveTime += dummyTime/1000.0; 
  nSolve++;
  filePrint(stderr,"Average Solve Time: %3.2e \n",reducedSolveTime/double(nSolve));

  // zero out unused parts of rhs
  for(int i=0; i<startCol_; ++i) rhs[i] = 0;
  for(int i=startCol_+blockCols_; i<rhs.size(); ++i) rhs[i] = 0;
}

template <typename Scalar, template<typename> class GenVecType, class BaseSolver>
void
GenEiSparseGalerkinProjectionSolver<Scalar,GenVecType,BaseSolver>::solve(GenVecType<Scalar> &rhs, GenVecType<Scalar> &sol)
{
  Eigen::Map< Eigen::Matrix<Scalar, Eigen::Dynamic, 1> > b(rhs.data()+startCol_, blockCols_), x(sol.data()+startCol_, blockCols_);
  sol.zero();

  double dummyTime = -1.0*getTime();
  if(dualBlockCols_ > 0 && contact_) {
    Eigen::Matrix<Scalar,Eigen::Dynamic,Eigen::Dynamic> CE(0,0), CI = -reducedConstraintMatrix_.block(startDualCol_,startCol_,dualBlockCols_,blockCols_).transpose();
    Eigen::Matrix<Scalar,Eigen::Dynamic,1> g0 = -b, ce0(0,1), _x(blockCols_), Lambda(0,1), Mu(dualBlockCols_);
    solve_quadprog2(llt_, c1_, g0, CE, ce0, CI, reducedConstraintRhs_.segment(startDualCol_,dualBlockCols_), _x, &Lambda, &Mu, tol_);
    x = _x;
    reducedConstraintForce_.setZero();
    reducedConstraintForce_.segment(startCol_,blockCols_) = reducedConstraintMatrix_.block(startDualCol_,startCol_,dualBlockCols_,blockCols_).transpose()*Mu;
  }
  else if(selfadjoint_ && !Empirical) x = llt_.solve(b);
  else x = lu_.solve(b);

  dummyTime += getTime();
  reducedSolveTime += dummyTime/1000.0;
  nSolve++;
  filePrint(stderr,"Average Solve Time: %3.2e \n",reducedSolveTime/double(nSolve));
}

template <typename Scalar, template<typename> class GenVecType, class BaseSolver>
double
GenEiSparseGalerkinProjectionSolver<Scalar,GenVecType,BaseSolver>::getResidualNorm(const GenVecType<Scalar> &v)
{
  Eigen::Map< Eigen::Matrix<Scalar, Eigen::Dynamic, 1> > vj(v.data()+startCol_, blockCols_);
  return vj.norm();
}

template <typename Scalar, template<typename> class GenVecType, class BaseSolver>
double
GenEiSparseGalerkinProjectionSolver<Scalar,GenVecType,BaseSolver>::getFNormSq(GenVecType<Scalar> &v)
{
  Eigen::Map< Eigen::Matrix<Scalar, Eigen::Dynamic, 1> > vj(v.data()+startCol_, blockCols_);
  double dummy = vj.norm();
  return dummy*dummy;
}

} /* end namespace Rom */

template class Rom::GenEiSparseGalerkinProjectionSolver<double, GenDistrVector, GenParallelSolver<double> >;
template class Rom::GenEiSparseGalerkinProjectionSolver<std::complex<double>, GenDistrVector, GenParallelSolver<std::complex<double> > >;

#endif /* USE_EIGEN3 */

#endif /* ROM_EIGALERKINPROJECTIONSOLVER_C */
