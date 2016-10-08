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

NonnegativeMatrixFactorization::NonnegativeMatrixFactorization(int maxBasisDimension, int method) :
  rowCount_(0),
  colCount_(0),
  basisDimension_(maxBasisDimension),
  maxBasisDimension_(maxBasisDimension),
  numRandInit_(1),
  method_(method),
  matrixBuffer_(0),
  robBuffer_(0),
  maxIter_(100)
{}

void
NonnegativeMatrixFactorization::matrixSizeIs(int rows, int cols) {
  matrixBuffer_.sizeIs(rows * cols);

  rowCount_ = rows;
  colCount_ = cols;
}

void
NonnegativeMatrixFactorization::robSizeIs(int rows, int cols) {
  robBuffer_.sizeIs(rows * cols);
  maxBasisDimensionIs(cols);
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
NonnegativeMatrixFactorization::maxBasisDimensionIs(int maxBasisDimension) {
  maxBasisDimension_ = maxBasisDimension;
}

void
NonnegativeMatrixFactorization::basisDimensionIs(int basisDimension) {
  basisDimension_ = basisDimension;
}

void 
NonnegativeMatrixFactorization::numRandInitIs(int numRandInit) {
  numRandInit_ = numRandInit;
}

void
NonnegativeMatrixFactorization::solve(int basisDimensionRestart) {
#ifdef USE_EIGEN3
  Eigen::Map<Eigen::MatrixXd> X(matrixBuffer_.array(), rowCount_, colCount_);

  std::vector<int> cols;
  for(int i=0; i<colCount_; ++i) if(!(X.col(i).array() == 0).all()) cols.push_back(i);
  std::cerr << "X has " << X.cols() << " columns of which " << cols.size() << " are non-zero\n";

  std::vector<int> rows;
  for(int i=0; i<rowCount_; ++i) if(!(X.row(i).array() == 0).all()) rows.push_back(i);
  std::cerr << "X has " << X.rows() << " rows of which " << rows.size() << " are non-zero\n";

  Eigen::Map<Eigen::MatrixXd> ROB(robBuffer_.array(), rowCount_, maxBasisDimension_);
  int m = rows.size();
  int n = cols.size();
  const int &k = basisDimension_;
  Eigen::MatrixXd A(m,n);
  for(int i=0; i<m; ++i){
    for(int j=0; j<n; ++j){
      A(i,j) = -X(rows[i],cols[j]); // note -ve is due to sign convention (Lagrange multipliers are negative in Aero-S)
    }
  }

  

  Eigen::MatrixXd W(m,k), H(k,n);

  double res, index;
  Eigen::MatrixXd Err(A);
  switch(method_) { 
    default: case 1 : { // NMF ROB
      Eigen::MatrixXd W_min(W);
      double res_min;
      for (int iRandom=0; iRandom<numRandInit_; ++iRandom) {
        // Initialization
        if (basisDimensionRestart>0) { 
          std::cout << "Initial guess based on smaller factorization" << std::endl;
          for (int i=0; i<m; ++i)
            for (int j=0; j<basisDimensionRestart; ++j)
              W(i,j) = ROB(rows[i],j);
          W.rightCols(k-basisDimensionRestart) = 0.5*(Eigen::MatrixXd::Random(m,k-basisDimensionRestart)+Eigen::MatrixXd::Ones(m,k-basisDimensionRestart));
        } else {
          W = 0.5*(Eigen::MatrixXd::Random(m,k) + Eigen::MatrixXd::Ones(m,k));
        }

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
        if (numRandInit_>1) {
          if (iRandom==0) {
            res_min = res;
            W_min = W;
          }
          else {
            if (res_min>res) {
              res_min = res;
              W_min = W;
            }
            if (iRandom==numRandInit_-1)
              W = W_min;
              std::cout << "NNMF basis retained with rel. residual = " << res_min << std::endl;
          }
        }
      }
    
    } break;
    case 2 : { // Greedy ROB
      W.setZero();
      H.setZero();
      // first vector is vector with largest magnitude
      index = findColumnWithLargestMagnitude(A);
      for (int i=1; i<=k; ++i) {
        W.col(i-1) = A.col(index);
        H.topRows(i) = solveNNLS_MRHS(W.leftCols(i), A);
        Err = A-W.leftCols(i)*H.topRows(i); 
        res = Err.norm()/A.norm();
        index = findColumnWithLargestMagnitude(Err);
        std::cout << "greedy iteration = " << i << ", rel. residual = " << res << ", maximum error = " << Err.col(index).norm() << std::endl;
      }
      // normalize vectors in W
      for (int i=0; i<k; ++i) W.col(i).normalize();
    } break;
  }

  // copy W into buffer
  ROB.setZero();
  for(int i=0; i<m; ++i)
    for(int j=0; j<k; ++j)
      ROB(rows[i],j) = W(i,j);
#endif
}

#ifdef USE_EIGEN3
Eigen::MatrixXd
NonnegativeMatrixFactorization::solveNNLS_MRHS(const Eigen::Ref<const Eigen::MatrixXd> &A, const Eigen::Ref<const Eigen::MatrixXd> &B)
{
  Eigen::MatrixXd X(A.cols(),B.cols());  
  for(int j=0; j<B.cols(); ++j) {
    X.col(j) = solveNNLS(A, B.col(j));
  }
  return X;
}

int
NonnegativeMatrixFactorization::findColumnWithLargestMagnitude(const Eigen::Ref<const Eigen::MatrixXd> &X)
{
  int index = 0;
  X.colwise().norm().maxCoeff(&index);
  return index;
}

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
