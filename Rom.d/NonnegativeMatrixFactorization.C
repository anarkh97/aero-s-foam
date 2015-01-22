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

  int NNMFROB = 0;
  double res, index;
  Eigen::MatrixXd Err(A);
  if (NNMFROB) { 
    W = Eigen::MatrixXd::Random(W.rows(), W.cols());
    Eigen::MatrixXd W_copy(W);
    for(int i=0; i<maxIter_; ++i) {
    
      H = solveNNLS_MRHS(W, A);
      W.transpose() = solveNNLS_MRHS(H.transpose(),A.transpose());
      Err = A-W*H;
      index = findColumnWithLargestMagnitude(Err);
      res = Err.norm()/A.norm();
      double inc = (W-W_copy).norm();
      std::cout << "iteration = " << i+1 << ", rel. residual = " << res << ", maximum error = " << Err.col(index).norm() << ", solution incr. = " << inc << std::endl;
      if(inc < tol_) break;
      W_copy = W;
    }
  }
  else {// Greedy ROB
    W = Eigen::MatrixXd::Zero(W.rows(), W.cols());
    H = Eigen::MatrixXd::Zero(H.rows(), H.cols());
    // first vector is vector with largest magnitude
    index = findColumnWithLargestMagnitude(A);
    W.col(0) = A.col(index);
    for (int i=1; i<k; ++i) {
      H.topRows(i) = solveNNLS_MRHS(W.leftCols(i), A);
      Err = A-W.leftCols(i)*H.topRows(i); 
      res = Err.norm()/A.norm();
      index = findColumnWithLargestMagnitude(Err);
      std::cout << "greedy iteration = " << i << ", rel. residual = " << res << ", maximum error = " << Err.col(index).norm() << std::endl;
      W.col(i) = A.col(index);
    }
    H.topRows(k) = solveNNLS_MRHS(W.leftCols(k), A);
    Err = A-W.leftCols(k)*H.topRows(k);
    res = Err.norm()/A.norm();
    index = findColumnWithLargestMagnitude(Err);
    std::cout << "greedy iteration = " << k << ", rel. residual = " << res << ", maximum error = " << Err.col(index).norm() << std::endl;
    //normalize vectors in W
    for (int i=0; i<k; ++i)
      W.col(i) = W.col(i)/W.col(i).norm();
  }

  // copy W into matrixBuffer
  X.setZero();
  for(int i=0; i<m; ++i)
    for(int j=0; j<k; ++j)
      X(rows[i],j) = W(i,j);
#endif
}

#ifdef USE_EIGEN3
Eigen::MatrixXd
NonnegativeMatrixFactorization::solveNNLS_MRHS(const Eigen::Ref<const Eigen::MatrixXd> &_A, const Eigen::Ref<const Eigen::MatrixXd> &_B)
{

  Eigen::MatrixXd X(_A.cols(),_B.cols());  
  for(int j=0; j<_B.cols(); ++j) {
    X.col(j) = solveNNLS(_A, _B.col(j));
   }
   return X;
}
#endif

#ifdef USE_EIGEN3
int
NonnegativeMatrixFactorization::findColumnWithLargestMagnitude(const Eigen::Ref<const Eigen::MatrixXd> &_X)
{

  int index = 0;
  double maxNorm = _X.col(0).norm();
  double vecNorm;
  for(int j=1; j<_X.cols(); ++j) {
    vecNorm = _X.col(j).norm();
    if (maxNorm < vecNorm) {
      index = j;
      maxNorm = vecNorm;
    }
  }
  return index;

}
#endif

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
