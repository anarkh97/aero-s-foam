#include "ParallelSparseNonNegativeLeastSquaresSolver.h"
#include <Timers.d/GetTime.h>
#include <Utils.d/DistHelper.h>

#include <stdexcept>
#include <vector>
#include <iostream>

#ifdef USE_EIGEN3
#include <Eigen/Core>

Eigen::Array<Eigen::VectorXd,Eigen::Dynamic,1>
pnncgp(const std::vector<Eigen::Map<Eigen::MatrixXd> >&A, const Eigen::Ref<const Eigen::VectorXd> &b, double& rnorm, const long int n,
       long int &info, double maxsze, double maxite, double reltol, bool verbose, bool scaling, double &dtime);

Eigen::Array<Eigen::VectorXd,Eigen::Dynamic,1>
pgpfp(const std::vector<Eigen::Map<Eigen::MatrixXd> >&A, const Eigen::Ref<const Eigen::VectorXd> &b, double& rnorm, const long int n,
      long int &info, double maxsze, double maxite, double reltol, bool verbose, bool scaling, bool positive, double &dtime);

Eigen::Array<Eigen::VectorXd,Eigen::Dynamic,1>
pomp(const std::vector<Eigen::Map<Eigen::MatrixXd> >&A, const Eigen::Ref<const Eigen::VectorXd> &b, double& rnorm, const long int n,
       long int &info, double maxsze, double maxite, double reltol, bool verbose, bool scaling, bool positive, double &dtime);

Eigen::Array<Eigen::VectorXd,Eigen::Dynamic,1>
splh(const std::vector<Eigen::Map<Eigen::MatrixXd> >&A, const Eigen::Ref<const Eigen::VectorXd> &b, double& rnorm, const long int n,
       long int &info, double maxsze, double maxite, double reltol, bool verbose, bool scaling, bool positive, double &dtime,
       int npMax, int scpkMB, int scpkNB, int scpkMP, int scpkNP);

#endif

namespace Rom {

ParallelSparseNonNegativeLeastSquaresSolver::ParallelSparseNonNegativeLeastSquaresSolver(int nsub, SparseNonNegativeLeastSquaresSolver<std::vector<double>,size_t> **sd) :
  equationCount_(0),
  unknownCount_(0),
  relativeTolerance_(1.0e-6),
  rhsBuffer_(0),
  errorMagnitude_(),
  verboseFlag_(true),
  scalingFlag_(true),
  projectFlag_(false),
  positivity_(true),
  solverType_(0),
  maxSizeRatio_(1.0),
  maxIterRatio_(3.0),
  npMax_(0),
  scpkMB_(0),
  scpkNB_(0),
  scpkMP_(0),
  scpkNP_(0),
  nsub_(nsub),
  sd_(sd)
{}

void
ParallelSparseNonNegativeLeastSquaresSolver::problemSizeIs(long eqnCount, long unkCount) {
  if (eqnCount < 0 || unkCount < 0) {
    throw std::domain_error("Illegal problem size");
  }

  equationCount_ = eqnCount;
  unknownCount_ = unkCount;
  rhsBuffer_.sizeIs(equationCount());
}

void
ParallelSparseNonNegativeLeastSquaresSolver::solve() {
  if (equationCount_ == 0 || unknownCount_ == 0) {
    return;
  }

  long int info;
  double dtime = 0.0;
  double t0 = getTime();
  switch(solverType_) {
    default :
    case 1 : { // Non-negative Conjugate Gradient Pursuit
#ifdef USE_EIGEN3
      filePrint(stderr, " ... Using Parallel NNCGP Solver    ...\n");
      std::vector<Eigen::Map<Eigen::MatrixXd> > A(nsub_, Eigen::Map<Eigen::MatrixXd>(NULL,0,0));
      Eigen::Array<Eigen::VectorXd,Eigen::Dynamic,1> x(nsub_);
      Eigen::Map<Eigen::VectorXd> b(rhsBuffer_.array(),equationCount_);
      for(int i=0; i<nsub_; ++i) {
        new (&A[i]) Eigen::Map<Eigen::MatrixXd>(&*sd_[i]->matrixBuffer(),sd_[i]->equationCount(),sd_[i]->unknownCount());
      }
      x = pnncgp(A, b, errorMagnitude_, unknownCount_, info, maxSizeRatio_, maxIterRatio_, relativeTolerance_, verboseFlag_, scalingFlag_, dtime);
      for(int i=0; i<nsub_; ++i) {
        Eigen::Map<Eigen::VectorXd>(const_cast<double*>(sd_[i]->solutionBuffer()),sd_[i]->unknownCount()) = x[i];
        A[i].~Map<Eigen::MatrixXd>();
      }
#else
      std::cerr << "USE_EIGEN3 is not defined here in ParallelSparseNonNegativeLeastSquaresSolver::solve\n";
      exit(-1);
#endif
    } break;

     case 2 : { // Polytope Faces Pursuit
#ifdef USE_EIGEN3
      filePrint(stderr, " ... Using Parallel GPFP Solver     ...\n");
      std::vector<Eigen::Map<Eigen::MatrixXd> > A(nsub_, Eigen::Map<Eigen::MatrixXd>(NULL,0,0));
      Eigen::Array<Eigen::VectorXd,Eigen::Dynamic,1> x(nsub_);
      Eigen::Map<Eigen::VectorXd> b(rhsBuffer_.array(),equationCount_);
      for(int i=0; i<nsub_; ++i) {
        new (&A[i]) Eigen::Map<Eigen::MatrixXd>(&*sd_[i]->matrixBuffer(),sd_[i]->equationCount(),sd_[i]->unknownCount());
      }
      x = pgpfp(A, b, errorMagnitude_, unknownCount_, info, maxSizeRatio_, maxIterRatio_, relativeTolerance_, verboseFlag_, scalingFlag_, positivity_, dtime);
      for(int i=0; i<nsub_; ++i) {
        Eigen::Map<Eigen::VectorXd>(const_cast<double*>(sd_[i]->solutionBuffer()),sd_[i]->unknownCount()) = x[i];
        A[i].~Map<Eigen::MatrixXd>();
      }
#else
      std::cerr << "USE_EIGEN3 is not defined here in ParallelSparseNonNegativeLeastSquaresSolver::solve\n";
      exit(-1);
#endif
    } break;
    
    case 5 : { // Orthogonal Matching Pursuit with Non-Negative L2 minimization
#ifdef USE_EIGEN3
      filePrint(stderr, " ... Using Parallel OMP Solver    ...\n");
      std::vector<Eigen::Map<Eigen::MatrixXd> > A(nsub_, Eigen::Map<Eigen::MatrixXd>(NULL,0,0));
      Eigen::Array<Eigen::VectorXd,Eigen::Dynamic,1> x(nsub_);
      Eigen::Map<Eigen::VectorXd> b(rhsBuffer_.array(),equationCount_);
      for(int i=0; i<nsub_; ++i) {
        new (&A[i]) Eigen::Map<Eigen::MatrixXd>(&*sd_[i]->matrixBuffer(),sd_[i]->equationCount(),sd_[i]->unknownCount());
      }
      x = pomp(A, b, errorMagnitude_, unknownCount_, info, maxSizeRatio_, maxIterRatio_, relativeTolerance_, verboseFlag_, scalingFlag_, positivity_, dtime);
      for(int i=0; i<nsub_; ++i) {
        Eigen::Map<Eigen::VectorXd>(const_cast<double*>(sd_[i]->solutionBuffer()),sd_[i]->unknownCount()) = x[i];
        A[i].~Map<Eigen::MatrixXd>();
      }
#else
      std::cerr << "USE_EIGEN3 is not defined here in ParallelSparseNonNegativeLeastSquaresSolver::solve\n";
      exit(-1);
#endif
    } break;

    case 6 : { // Scalapack LH solver
#ifdef USE_EIGEN3
      filePrint(stderr, " ... Using Scalapack LH Solver    ...\n");
      std::vector<Eigen::Map<Eigen::MatrixXd> > A(nsub_, Eigen::Map<Eigen::MatrixXd>(NULL,0,0));
      Eigen::Array<Eigen::VectorXd,Eigen::Dynamic,1> x(nsub_);
      Eigen::Map<Eigen::VectorXd> b(rhsBuffer_.array(),equationCount_);
      for(int i=0; i<nsub_; ++i) {
        new (&A[i]) Eigen::Map<Eigen::MatrixXd>(&*sd_[i]->matrixBuffer(),sd_[i]->equationCount(),sd_[i]->unknownCount());
      }
      x = splh(A, b, errorMagnitude_, unknownCount_, info, maxSizeRatio_, maxIterRatio_, relativeTolerance_, verboseFlag_, scalingFlag_, projectFlag_, dtime, npMax_, scpkMB_, scpkNB_, scpkMP_, scpkNP_);
      for(int i=0; i<nsub_; ++i) {
        Eigen::Map<Eigen::VectorXd>(const_cast<double*>(sd_[i]->solutionBuffer()),sd_[i]->unknownCount()) = x[i];
        A[i].~Map<Eigen::MatrixXd>();
      }
#else
      std::cerr << "USE_EIGEN3 is not defined here in ParallelSparseNonNegativeLeastSquaresSolver::solve\n";
      exit(-1);
#endif
    }
  }
  double t = (getTime() - t0)/1000.0;
  filePrint(stderr, " ... Solve Time = %13.6f s   ...\n",t);
  filePrint(stderr, " ... DDate Time = %13.6f s   ...\n",dtime);
  filePrint(stderr, " ... %% DDate   =  %13.6f %% ...\n",(dtime/t)*100.0);

  if (info == 2) {
    throw std::logic_error("Illegal problem size");
  }

  if (info == 3) {
    throw std::runtime_error("Solution did not converge");
  }
}

} // end namespace Rom
