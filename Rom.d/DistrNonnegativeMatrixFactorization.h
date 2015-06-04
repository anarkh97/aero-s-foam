#ifndef ROM_DISTRNONNEGATIVEMATRIXFACTORIZATION_H
#define ROM_DISTRNONNEGATIVEMATRIXFACTORIZATION_H

class Communicator;
class SCDoubleMatrix;

#include <Eigen/Core>

namespace Rom {

class DistrNonnegativeMatrixFactorization {
public:
  // Local data distribution
  int localRows() const { return localRows_; }

  // Local buffers: Internal column-major ordering, zero-based indexing
  // Local matrix buffer: [localRows by colCount]
  double *matrixColBuffer(int col);
  // Local basis buffer: [localRows by basisDimension]
  const double *basisColBuffer(int col) const;

  void solve(int);

  DistrNonnegativeMatrixFactorization(Communicator * comm, int rowCount, int colCount, int localRows, int basisDimension, int blockSize, int maxIter, double tol);
  ~DistrNonnegativeMatrixFactorization();

private:
  void solveNNLS_MRHS(SCDoubleMatrix &A, SCDoubleMatrix &B, SCDoubleMatrix &X, int flag);

  // Disallow copy & assignment
  DistrNonnegativeMatrixFactorization(const DistrNonnegativeMatrixFactorization &);
  DistrNonnegativeMatrixFactorization & operator=(const DistrNonnegativeMatrixFactorization &);

  Communicator * communicator_;
  int rowCount_, colCount_, localRows_, basisDimension_, blockSize_, maxIter_;
  double tol_;

  Eigen::MatrixXd matrixBuffer_;
  Eigen::MatrixXd basisBuffer_;
};

/* Helper functions for buffer access */
inline
double *
DistrNonnegativeMatrixFactorization::matrixColBuffer(int col) {
  return matrixBuffer_.data() + col*localRows_;
}

inline
const double *
DistrNonnegativeMatrixFactorization::basisColBuffer(int col) const {
  return basisBuffer_.data() + col*localRows_;
}

} // end namespace Rom

#endif /* ROM_DISTRNONNEGATIVEMATRIXFACTORIZATION_H */
