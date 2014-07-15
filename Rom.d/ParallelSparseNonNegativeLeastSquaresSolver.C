#include "ParallelSparseNonNegativeLeastSquaresSolver.h"
#include <Utils.d/DistHelper.h>

#include <stdexcept>

#ifdef USE_EIGEN3
#include <Eigen/Core>

Eigen::Array<Eigen::VectorXd,Eigen::Dynamic,1>
pnncgp(const Eigen::Array<Eigen::MatrixXd,Eigen::Dynamic,1> &A, const Eigen::Ref<const Eigen::VectorXd> &b, double& rnorm, const long int n,
       double maxsze, double reltol, bool verbose, bool scaling);
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
  solverType_(0),
  maxSizeRatio_(1.0),
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

  switch(solverType_) {
    default :
    case 1 : { // Non-negative Conjugate Gradient Pursuit
#ifdef USE_EIGEN3
      filePrint(stderr, " ... Using Parallel NNCGP Solver    ...\n");
      Eigen::Array<Eigen::MatrixXd,Eigen::Dynamic,1> A(nsub_);
      Eigen::Array<Eigen::VectorXd,Eigen::Dynamic,1> x(nsub_);
      Eigen::Map<Eigen::VectorXd> b(rhsBuffer_.array(),equationCount_);
      for(int i=0; i<nsub_; ++i) {
        A[i] = Eigen::Map<Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,Eigen::ColMajor> >(&*sd_[i]->matrixBuffer(),sd_[i]->equationCount(),sd_[i]->unknownCount());
      }
      x = pnncgp(A, b, errorMagnitude_, unknownCount_, maxSizeRatio_, relativeTolerance_, verboseFlag_, scalingFlag_);
      for(int i=0; i<nsub_; ++i) {
        Eigen::Map<Eigen::VectorXd>(const_cast<double*>(sd_[i]->solutionBuffer()),sd_[i]->unknownCount()) = x[i];
      }
#else
      std::cerr << "USE_EIGEN3 is not defined here in ParallelSparseNonNegativeLeastSquaresSolver::solve\n";
      exit(-1);
#endif
    }
  }
}

} // end namespace Rom
