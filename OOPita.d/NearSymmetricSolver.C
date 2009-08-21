#include "NearSymmetricSolver.h"

#include <algorithm>
#include <iterator>

extern "C" {
  // Lapack: Matrix-matrix multiplication
  // C := alpha*op(A)*op(B) + beta*C where op(M) = M or M^T
  /*void _FORTRAN(dgemm)(const char* transa, const char* transb, const int* m, const int* n, const int* k,
                       const double* alpha, const double* a, const int* lda, const double* b, const int* ldb,
                       const double* beta, double* c, const int* ldc);*/

  // Blas: Rank-1 update A += alpha*x*y'
  void _FORTRAN(dger)(const int* m, const int* n, const double* alpha,
                      const double* x, const int* incx, const double* y, const int* incy,
                      double* a, const int* lda);

  // Blas: Triangular system resolution
  void _FORTRAN(dtrsv)(const char * uplo, const char * trans, const char * diag, const int * n,
                       const double * a, const int * lda, double * x, const int * incx);

  // Lapack: Row pivoting
  void _FORTRAN(dlaswp)(const int* n, double* a, const int* lda, const int* k1, const int* k2,
                        const int* ipiv, const int* incx);

  // Blas: Column (vector) interchange
  void _FORTRAN(dswap)(const int* n, double* x, const int* incx, double* y, const int* incy);

  // Blas: Vector scaling
  void _FORTRAN(dscal)(const int* n, const double* da, double* dx, const int* incx);

  // Lapack: Diagonal scaling
  void _FORTRAN(dpoequ)(const int* n, const double* a, const int* lda,
                        double* s, double* scond, double* amax, int* info);

  // Lapack: Perform scaling
  void _FORTRAN(dlaqge)(const int* m, const int* n, double* a, const int* lda,
                        const double* r, const double* c, const double* rowcnd, const double* colcnd,
                        const double* amax, char* equed);
}

namespace Pita {

NearSymmetricSolver::NearSymmetricSolver(double tol) :
  RankDeficientSolver(),
  tolerance_(tol),
  transposedMatrix_(),
  pivots_(),
  rescalingStatus_(NO_RESCALING),
  scaling_()
{}

void
NearSymmetricSolver::toleranceIs(double tol) {
  // TODO Extend/shrink factorization ?
  setTolerance(tol);
}

void
NearSymmetricSolver::transposedMatrixIs(const FullSquareMatrix & tm) {
  transposedMatrix_.copy(tm);
  
  setMatrixSize(transposedMatrix_.dim());
  setFactorRank(0);
  setStatus(NON_FACTORIZED);
  rescalingStatus_ = NO_RESCALING;
}

void
NearSymmetricSolver::statusIs(RankDeficientSolver::Status s) {
  if (status() == s) {
    return;
  }

  if (s == FACTORIZED) {
    // Rescale matrix
    scaling_.sizeIs(matrixSize());
    int info;
    double scond, amax;
    _FORTRAN(dpoequ)(&getMatrixSize(), transposedMatrix_.data(), &getMatrixSize(),
                    scaling_.array(), &scond, &amax, &info);
    char equed;
    _FORTRAN(dlaqge)(&getMatrixSize(), &getMatrixSize(), transposedMatrix_.data(), &getMatrixSize(),
                     scaling_.array(), scaling_.array(), &scond, &scond, &amax, &equed); 

    switch (equed) {
      case 'N':
        rescalingStatus_ = NO_RESCALING;
        log() << "No rescaling\n";
        break;
      case 'R':
        rescalingStatus_ = ROW_RESCALING;
        log() << "Row rescaling\n";
        break;
      case 'B':
        rescalingStatus_ = SYMMETRIC_RESCALING;
        log() << "Symmetric rescaling\n";
        break;
      default:
        throw Fwk::InternalException("NearSymmetricSolver::statusIs - Invalid rescaling status");
    }

    // Initialize permutation
    getFactorPermutation().sizeIs(matrixSize());
    for (int i = 0; i < matrixSize(); ++i) {
      getFactorPermutation()[i] = i + 1;
    }
    pivots_.sizeIs(matrixSize());

    // Initialize diagonal pivot values
    SimpleBuffer<double> pivot_values;
    pivot_values.sizeIs(matrixSize());

    // Numerical constant
    const int int_one = 1;
    const double minus_one = -1;

    double first_pivot = 0.0;

    int k;
    for (k = 0; k < matrixSize(); ++k) {
      // Update pivot values, find largest and permute
      for (int i = 0; i < matrixSize(); ++i) {
        pivot_values[i] = transposedMatrix()[i][i];
      }
      double * head_pivot = pivot_values.array() + k;
      double * max_pivot = std::max_element(head_pivot, pivot_values.array() + matrixSize());
      const int p = std::distance(pivot_values.array(), max_pivot);
      pivots_[k] = p + 1;
      std::swap(head_pivot, max_pivot);
      std::swap(getFactorPermutation()[k], getFactorPermutation()[p]);
      
      if (k == 0) {
        first_pivot = *head_pivot;
      }

      // Perform symmetric permutation
      _FORTRAN(dswap)(&getMatrixSize(), &transposedMatrix_[k][0], &int_one, &transposedMatrix_[p][0], &int_one);
      _FORTRAN(dswap)(&getMatrixSize(), &transposedMatrix_[0][k], &getMatrixSize(), &transposedMatrix_[0][p], &getMatrixSize());
      
      // Check for singularity
      if (*head_pivot <= first_pivot * tolerance()) {
        break;
      }

      // Compute column of L
      const double pivot_inverse = 1.0 / *head_pivot;
      const int remainder_size = matrixSize() - (k + 1);
      _FORTRAN(dscal)(&remainder_size, &pivot_inverse, &transposedMatrix_[k][k+1], &int_one);

      // Update remaining block of A
      _FORTRAN(dger)(&remainder_size, &remainder_size, &minus_one,
                     &transposedMatrix_[k][k+1], &int_one, &transposedMatrix_[k+1][k], &getMatrixSize(),
                     &transposedMatrix_[k+1][k+1], &getMatrixSize());
    }

    setFactorRank(k);
  }

  setStatus(s);
}

const Vector &
NearSymmetricSolver::solution(Vector & rhs) const {
  if (rhs.size() != matrixSize()) {
    throw Fwk::RangeException("in NearSymmetricSolver::solution - Size mismatch"); 
  }

  if (status() != FACTORIZED) {
    throw Fwk::RangeException("in NearSymmetricSolver::solution - Non-factorized matrix");
  }

  // 1) rhs <- P * rhs
  const int int_one = 1;
  _FORTRAN(dlaswp)(&int_one, rhs.data(), &getMatrixSize(), &int_one, &getFactorRank(),
                   pivots_.array(), &int_one);

  // 2) rhs <- S^{-1} * rhs
  if (rescalingStatus_ != NO_RESCALING) {
    for (int i = 0; i < factorRank(); ++i) {
      rhs[i] *= scaling_[factorPermutation(i)];
    }
  }

  // 3) rhs <- L^{-1} * rhs
  const char lower = 'L';
  const char non_trans = 'N';
  const char unit_diag = 'U';
  _FORTRAN(dtrsv)(&lower, &non_trans, &unit_diag, &getFactorRank(), transposedMatrix().data(),
                  &getMatrixSize(), rhs.data(), &int_one);

  // 4) rhs <- U^{-1} * rhs
  const char upper = 'U';
  const char non_unit_diag = 'N';
  _FORTRAN(dtrsv)(&upper, &non_trans, &non_unit_diag, &getFactorRank(), transposedMatrix().data(),
                  &getMatrixSize(), rhs.data(), &int_one);
  
  // 5) rhs <- S^{-1} * rhs
  if (rescalingStatus_ == SYMMETRIC_RESCALING) {
    for (int i = 0; i < factorRank(); ++i) {
      rhs[i] *= scaling_[factorPermutation(i)];
    }
  }

  // 6) Pad with zeros
  std::fill(rhs.data() + getFactorRank(), rhs.data() + getMatrixSize(), 0.0);

  // 7) rhs <- P * rhs
  const int int_minus_one = -1;
  _FORTRAN(dlaswp)(&int_one, rhs.data(), &getMatrixSize(), &int_one, &getFactorRank(),
                   pivots_.array(), &int_minus_one);

  return rhs;
}

} /* end namespace Pita */
