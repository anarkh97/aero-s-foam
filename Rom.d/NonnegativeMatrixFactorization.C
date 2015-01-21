#include "NonnegativeMatrixFactorization.h"

#include <Utils.d/linkfc.h>

#include <iostream>
#include <vector>

extern "C" {
  void _FORTRAN(spnnls)(double *a, const long int *mda, const long int *m, const long int *n,
                        const double *b, double *x, const double *reltol, double *rnorm, double *w,
                        double *zz, double *zz2, long int *index, long int *mode, long int *prtflg,
                        long int *sclflg, const double *maxsze, const double *maxite, double *dtime);
}

namespace Rom {

NonnegativeMatrixFactorization::NonnegativeMatrixFactorization(int basisDimension) :
  rowCount_(0),
  colCount_(0),
  basisDimension_(basisDimension),
  matrixBuffer_(0),
  maxIter_(100)
{}

void
NonnegativeMatrixFactorization::matrixSizeIs(int rows, int cols) {
  matrixBuffer_.sizeIs(rows * cols);

  rowCount_ = rows;
  colCount_ = cols;
}

void
NonnegativeMatrixFactorization::maxIterIs(int maxIter) {
  maxIter_ = maxIter;
}

void
NonnegativeMatrixFactorization::toleranceIs(double tol) {
  tol_ = tol;
}

void
NonnegativeMatrixFactorization::solve() {
#ifdef USE_EIGEN3
  Eigen::Map<Eigen::MatrixXd> X(matrixBuffer_.array(), rowCount_, colCount_);

  std::vector<int> cols;
  for(int i=0; i<colCount_; ++i) if(!(X.col(i).array() == 0).all()) cols.push_back(i);
  std::cerr << "X has " << X.cols() << " columns of which " << cols.size() << " are non-zero\n";

  std::vector<int> rows;
  for(int i=0; i<rowCount_; ++i) if(!(X.row(i).array() == 0).all()) rows.push_back(i);
  std::cerr << "X has " << X.rows() << " rows of which " << rows.size() << " are non-zero\n";

  int m = rows.size();
  int n = cols.size();
  const int &k = basisDimension_;
  Eigen::MatrixXd A(m,n);
  for(int i=0; i<m; ++i)
    for(int j=0; j<n; ++j)
      A(i,j) = -X(rows[i],cols[j]); // note -ve is due to sign convention (Lagrange multipliers are negative in Aero-S)
  Eigen::MatrixXd W(m,k), H(k,n);

  W = Eigen::MatrixXd::Random(W.rows(), W.cols());
  Eigen::MatrixXd W_copy(W);
  for(int i=0; i<maxIter_; ++i) {
    for(int j=0; j<A.cols(); ++j) {
      H.col(j) = solveNNLS(W, A.col(j));
    }
    for(int j=0; j<A.rows(); ++j) {
      W.row(j).transpose() = solveNNLS(H.transpose(), A.row(j).transpose());
    }
    double res = (A-W*H).norm()/A.norm();
    double inc = (W-W_copy).norm();
    std::cout << "iteration = " << i+1 << ", rel. residual = " << res << ", solution incr. = " << inc << std::endl;
    if(inc < tol_) break;
    W_copy = W;
  }

  // copy W into matrixBuffer
  X.setZero();
  for(int i=0; i<m; ++i)
    for(int j=0; j<k; ++j)
      X(rows[i],j) = W(i,j);
#endif
}

#ifdef USE_EIGEN3
Eigen::VectorXd
NonnegativeMatrixFactorization::solveNNLS(const Eigen::Ref<const Eigen::MatrixXd> &_A, const Eigen::Ref<const Eigen::VectorXd> &_b)
{
  Eigen::MatrixXd A(_A);
  const long int m = A.rows();
  const long int n = A.cols();
  Eigen::VectorXd b(_b);
  Eigen::VectorXd x(n);
  const double reltol = 1e-16;
  double rnorm;
  Eigen::Array<double,Eigen::Dynamic,1> w(n);
  Eigen::Array<double,Eigen::Dynamic,1> zz(m);
  Eigen::Array<double,Eigen::Dynamic,1> zz2(n);
  Eigen::Array<long int,Eigen::Dynamic,1> index(n);
  long int mode;
  long int prtflg = 0;
  long int scaflg = 1;
  const double maxsze = 1.0;
  const double maxite = 3.0;
  double dtime;

  _FORTRAN(spnnls)(A.data(), &m, &m, &n, b.data(), x.data(), &reltol, &rnorm, w.data(),
                   zz.data(), zz2.data(), index.data(), &mode, &prtflg, &scaflg, &maxsze, &maxite, &dtime);

  if(mode != 1) std::cerr << "Error: spnnls unsuccessful, mode = " << mode << std::endl;

  return x;
}
#endif

} /* end namespace Rom */
