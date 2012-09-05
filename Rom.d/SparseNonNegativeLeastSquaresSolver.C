#include "SparseNonNegativeLeastSquaresSolver.h"

#include <Utils.d/linkfc.h>

#include <stdexcept>
#ifdef USE_STXXL
#include "stxxl_matrix2d.hpp"
#else
#include <vector>
#endif

/*
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
  void _FORTRAN(spnnls)(double *a, const int *mda, const int *m, const int *n,
                        double *b, double *x, const double *reltol, double *rnorm, double *w,
                        double *zz, double *zz2, int *index, int *mode);
}
*/
#include "LawsonHanson.d/spnnls.cpp"

namespace Rom {

SparseNonNegativeLeastSquaresSolver::SparseNonNegativeLeastSquaresSolver() :
  equationCount_(0),
  unknownCount_(0),
  matrixLeadDim_(0),
  relativeTolerance_(1.0e-6),
  matrixBuffer_(),
  rhsBuffer_(0),
  solutionBuffer_(0),
  dualSolutionBuffer_(0),
  errorMagnitude_()
{}

void
SparseNonNegativeLeastSquaresSolver::problemSizeIs(int eqnCount, int unkCount) {
  if (eqnCount < 0 || unkCount < 0) {
    throw std::domain_error("Illegal problem size");
  }

  equationCount_ = matrixLeadDim_ = eqnCount;
  unknownCount_ = unkCount;
  matrixBuffer_.resize(matrixLeadDim_ * unknownCount());
  rhsBuffer_.sizeIs(equationCount());
  solutionBuffer_.sizeIs(unknownCount());
  dualSolutionBuffer_.sizeIs(unknownCount());
}

void
SparseNonNegativeLeastSquaresSolver::solve() {
  if (matrixBuffer_.size() == 0) {
    return;
  }

  SimpleBuffer<Scalar> workspace(equationCount());
  SimpleBuffer<Scalar> workspace2(unknownCount());
  SimpleBuffer<int> index(unknownCount());
  int info;

  stxxl_matrix2d<MatrixBufferType> A(&matrixBuffer_, matrixLeadDim_, unknownCount_);

  spnnls(A, matrixLeadDim_, equationCount_, unknownCount_,
         rhsBuffer_.array(), solutionBuffer_.array(), relativeTolerance_, errorMagnitude_, dualSolutionBuffer_.array(),
         workspace.array(), workspace2.array(), index.array(), info);

  if (info == 2) {
    throw std::logic_error("Illegal problem size");
  }

  if (info == 3) {
    throw std::runtime_error("Solution did not converge");
  }
}

} // end namespace Rom
