#include "SparseNonNegativeLeastSquaresSolver.h"
#include <Timers.d/GetTime.h>
#include <Utils.d/linkfc.h>

#include <stdexcept>
#include <cstdio>
#include <iostream>

extern "C" {
  // Approximately solve the sparse non-negative least-squares problem
  //   min support(x) st ||A * x - b|| < reltol * ||b|| and x >= 0
  // Input: A is (mda x n), b is (m x 1), reltol is scalar
  // Output: A <- Q A, b <- Q b where Q is (m x m) orthogonal,
  //         rnorm <- ||b - Ax||_2,
  //         x is the (n x 1) primal solution,
  //         w <- A^T(b - Ax) is the (n x 1) dual solution
  // Work: zz is (m x 1), zz2 is (n x 1), index is (n x 1)
  // Info: mode: 1 => success, 2 => bad dim, 3 => too many iter
  void _FORTRAN(spnnls)(double *a, const long int *mda, const long int *m, const long int *n,
                        double *b, double *x, const double *reltol, double *rnorm, double *w,
                        double *zz, double *zz2, long int *index, long int *mode, long int *prtflg,
                        long int *sclflg, const double *maxsze, const double *maxite, double *dtime);
}

#ifdef USE_EIGEN3
#include <Eigen/Core>

Eigen::VectorXd
nncgp(const Eigen::Ref<const Eigen::MatrixXd> &A, const Eigen::Ref<const Eigen::VectorXd> &b, double& rnorm,
      long int &info, double maxsze, double maxite, double reltol, bool verbose, bool scaling, double &dtime);

Eigen::VectorXd
gpfp(const Eigen::Ref<const Eigen::MatrixXd> &A, const Eigen::Ref<const Eigen::VectorXd> &b, double& rnorm,
     long int &info, double maxsze, double maxite, double reltol, bool verbose, bool scaling, bool positive, double &dtime);

Eigen::VectorXd
lars(const Eigen::Ref<const Eigen::MatrixXd> &A, const Eigen::Ref<const Eigen::VectorXd> &b, double& rnorm,
     long int &info, double maxsze, double maxite, double reltol, bool verbose, bool scaling, bool project, bool positive, double &dtime);

Eigen::VectorXd
mp(const Eigen::Ref<const Eigen::MatrixXd> &A, const Eigen::Ref<const Eigen::VectorXd> &b, double& rnorm,
     long int &info, double maxsze, double maxite, double reltol, bool verbose, bool scaling, bool positive);

Eigen::VectorXd
omp(const Eigen::Ref<const Eigen::MatrixXd> &A, const Eigen::Ref<const Eigen::VectorXd> &b, double& rnorm,
      long int &info, double maxsze, double maxite, double reltol, bool verbose, bool scaling, bool positive, double &dtime);
#endif

namespace Rom {

template<typename MatrixBufferType, typename SizeType>
SparseNonNegativeLeastSquaresSolver<MatrixBufferType,SizeType>::SparseNonNegativeLeastSquaresSolver() :
  equationCount_(0),
  unknownCount_(0),
  relativeTolerance_(1.0e-6),
  matrixBuffer_(),
  rhsBuffer_(0),
  solutionBuffer_(0),
  dualSolutionBuffer_(0),
  errorMagnitude_(),
  verboseFlag_(true),
  scalingFlag_(true),
  projectFlag_(false),
  positivity_(true),
  solverType_(0),
  maxSizeRatio_(1.0),
  maxIterRatio_(3.0)
{}

template<typename MatrixBufferType, typename SizeType>
void
SparseNonNegativeLeastSquaresSolver<MatrixBufferType,SizeType>::problemSizeIs(long eqnCount, long unkCount) {
  if (eqnCount < 0 || unkCount < 0) {
    throw std::domain_error("Illegal problem size");
  }

  SizeType bufSize = (eqnCount) * (unkCount);

  equationCount_ = eqnCount;
  unknownCount_ = unkCount;
  matrixBuffer_.resize(bufSize);
  rhsBuffer_.sizeIs(equationCount());
  solutionBuffer_.sizeIs(unknownCount());
  dualSolutionBuffer_.sizeIs(unknownCount());
}

template<typename MatrixBufferType, typename SizeType>
void
SparseNonNegativeLeastSquaresSolver<MatrixBufferType,SizeType>::solve() {
  if (matrixBuffer_.size() == 0) {
    return;
  }

  long int info;
  double dtime = 0.; 
  double t0 = getTime();
  switch(solverType_) {
    default :
    case 0 : { // Lawson & Hanson
      fprintf(stderr, " ... Using Lawson & Hanson Solver   ...\n");
      SimpleBuffer<Scalar> workspace(equationCount());
      SimpleBuffer<Scalar> workspace2(unknownCount());
      SimpleBuffer<long int> index(unknownCount());
      long int prtflg = (verboseFlag_) ? 1 : 0;
      long int scaflg = (scalingFlag_) ? 1 : 0;

      _FORTRAN(spnnls)(matrixBuffer_.data(), &equationCount_, &equationCount_, &unknownCount_, rhsBuffer_.array(),
                       solutionBuffer_.array(), &relativeTolerance_, &errorMagnitude_, dualSolutionBuffer_.array(),
                       workspace.array(), workspace2.array(), index.array(), &info, &prtflg, &scaflg, &maxSizeRatio_, &maxIterRatio_, &dtime);
    } break;

    case 1 : { // Non-negative Conjugate Gradient Pursuit
#ifdef USE_EIGEN3
      fprintf(stderr, " ... Using NNCGP Solver             ...\n");
      Eigen::Map<Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,Eigen::ColMajor> > A(matrixBuffer_.data(),equationCount_,unknownCount_);
      Eigen::Map<Eigen::VectorXd> x(solutionBuffer_.array(),unknownCount_);
      Eigen::Map<Eigen::VectorXd> b(rhsBuffer_.array(),equationCount_);
      x = nncgp(A, b, errorMagnitude_, info, maxSizeRatio_, maxIterRatio_, relativeTolerance_, verboseFlag_, scalingFlag_, dtime);
#else
      std::cerr << "USE_EIGEN3 is not defined here in SparseNonNegativeLeastSquaresSolver::solve\n";
      exit(-1);
#endif
    } break;
 
    case 2 : { // Gradient Polytope Faces Pursuit Pursuit
#ifdef USE_EIGEN3
      fprintf(stderr, " ... Using GPFP Solver              ...\n");
      Eigen::Map<Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,Eigen::ColMajor> > A(matrixBuffer_.data(),equationCount_,unknownCount_);
      Eigen::Map<Eigen::VectorXd> x(solutionBuffer_.array(),unknownCount_);
      Eigen::Map<Eigen::VectorXd> b(rhsBuffer_.array(),equationCount_);
      x = gpfp(A, b, errorMagnitude_, info, maxSizeRatio_, maxIterRatio_, relativeTolerance_, verboseFlag_, scalingFlag_, positivity_, dtime);
#else
      std::cerr << "USE_EIGEN3 is not defined here in SparseNonNegativeLeastSquaresSolver::solve\n";
      exit(-1);
#endif
    } break;

    case 3 : { // Least Angle Regression with LASSO modification
#ifdef USE_EIGEN3
      fprintf(stderr, " ... Using LARS Solver              ...\n");
      Eigen::Map<Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,Eigen::ColMajor> > A(matrixBuffer_.data(),equationCount_,unknownCount_);
      Eigen::Map<Eigen::VectorXd> x(solutionBuffer_.array(),unknownCount_);
      Eigen::Map<Eigen::VectorXd> b(rhsBuffer_.array(),equationCount_);
      x = lars(A, b, errorMagnitude_, info, maxSizeRatio_, maxIterRatio_, relativeTolerance_, verboseFlag_, scalingFlag_, projectFlag_, positivity_, dtime);
#else
      std::cerr << "USE_EIGEN3 is not defined here in SparseNonNegativeLeastSquaresSolver::solve\n";
      exit(-1);
#endif
    } break;
    case 4 : { // Matching Pursuits
#ifdef USE_EIGEN3
      fprintf(stderr, " ... Using Matching Pursuits Solver              ...\n");
      Eigen::Map<Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,Eigen::ColMajor> > A(matrixBuffer_.data(),equationCount_,unknownCount_);
      Eigen::Map<Eigen::VectorXd> x(solutionBuffer_.array(),unknownCount_);
      Eigen::Map<Eigen::VectorXd> b(rhsBuffer_.array(),equationCount_);
      x = mp(A, b, errorMagnitude_, info, maxSizeRatio_, maxIterRatio_, relativeTolerance_, verboseFlag_, scalingFlag_, positivity_);
#else
      std::cerr << "USE_EIGEN3 is not defined here in SparseNonNegativeLeastSquaresSolver::solve\n";
      exit(-1);
#endif
    } break;
    case 5 : { // Orthogonal Matching Pursuits
#ifdef USE_EIGEN3
      fprintf(stderr, " ... Using Orthogonal Matching Pursuits Solver              ...\n");
      Eigen::Map<Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,Eigen::ColMajor> > A(matrixBuffer_.data(),equationCount_,unknownCount_);
      Eigen::Map<Eigen::VectorXd> x(solutionBuffer_.array(),unknownCount_);
      Eigen::Map<Eigen::VectorXd> b(rhsBuffer_.array(),equationCount_);
      x = omp(A, b, errorMagnitude_, info, maxSizeRatio_, maxIterRatio_, relativeTolerance_, verboseFlag_, scalingFlag_, positivity_, dtime);
#else
      std::cerr << "USE_EIGEN3 is not defined here in SparseNonNegativeLeastSquaresSolver::solve\n";
      exit(-1);
#endif
    } break;

  }
  double t = (getTime() - t0)/1000.0;
  fprintf(stderr, " ... Solve Time    = %13.6f s   ...\n",t);
  fprintf(stderr, " ... DDate Time    = %13.6f s   ...\n",dtime);
  fprintf(stderr, " ... %% DDate       = %13.6f %% ...\n",(dtime/t)*100.);
  if (info == 2) {
    throw std::logic_error("Illegal problem size");
  }

  if (info == 3) {
    throw std::runtime_error("Solution did not converge");
  }
}

} // end namespace Rom
