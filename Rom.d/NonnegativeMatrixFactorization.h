#ifndef ROM_NONNEGATIVEMATRIXFACTORIZATION_H
#define ROM_NONNEGATIVEMATRIXFACTORIZATION_H

#include "SimpleBuffer.h"

#ifdef USE_EIGEN3
#include <Eigen/Core>
#endif

namespace Rom {

class NonnegativeMatrixFactorization {
public:
  int rowCount()           const { return rowCount_; }
  int colCount()           const { return colCount_; }

  void matrixSizeIs(int rows, int cols);
  void maxIterIs(int maxIter);
  void toleranceIs(double tol);

  const double *matrixBuffer() const;
  const double *matrixCol(int c) const;
  const double &matrixEntry(int r, int c) const;
  
  double *matrixBuffer();
  double *matrixCol(int c);
  double &matrixEntry(int r, int c);

  void solve();

  NonnegativeMatrixFactorization(int basisDimension, int method);

private:
  int rowCount_, colCount_, basisDimension_, method_, maxIter_;
  double tol_;

  SimpleBuffer<double> matrixBuffer_;

  // Disallow copy and assignment
  NonnegativeMatrixFactorization(const NonnegativeMatrixFactorization &);
  NonnegativeMatrixFactorization &operator=(const NonnegativeMatrixFactorization &);

#ifdef USE_EIGEN3
  Eigen::VectorXd solveNNLS(const Eigen::Ref<const Eigen::MatrixXd> &_A, const Eigen::Ref<const Eigen::VectorXd> &b);
  Eigen::MatrixXd solveNNLS_MRHS(const Eigen::Ref<const Eigen::MatrixXd> &_A, const Eigen::Ref<const Eigen::MatrixXd> &B);
  int findColumnWithLargestMagnitude(const Eigen::Ref<const Eigen::MatrixXd> &_X);
#endif
};

inline
const double *
NonnegativeMatrixFactorization::matrixBuffer() const {
  return matrixBuffer_.array();
}

inline
const double *
NonnegativeMatrixFactorization::matrixCol(int c) const {
  return matrixBuffer() + (c * rowCount_);
}

inline
const double &
NonnegativeMatrixFactorization::matrixEntry(int r, int c) const {
  return *(matrixCol(c) + r);
}

inline
double *
NonnegativeMatrixFactorization::matrixBuffer() {
  return const_cast<double *>(const_cast<const NonnegativeMatrixFactorization *>(this)->matrixBuffer());
}

inline
double *
NonnegativeMatrixFactorization::matrixCol(int c) {
  return const_cast<double *>(const_cast<const NonnegativeMatrixFactorization *>(this)->matrixCol(c));
}

inline
double &
NonnegativeMatrixFactorization::matrixEntry(int r, int c) {
  return const_cast<double &>(const_cast<const NonnegativeMatrixFactorization *>(this)->matrixEntry(r, c));
}

} /* end namespace Rom */

#endif /* ROM_NONNEGATIVEMATRIXFACTORIZATION_H */
