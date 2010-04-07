#include "PivotedCholeskySolver.h"
#include <algorithm>

/* Kernel routines */

extern "C" {
  // Pivoted Cholesky (LAPACK 3.2 routine calling BLAS level 3, also in lev3pchol.f)
  void _FORTRAN(dpstrf)(const char* uplo, const int* n, double* a, const int* lda, int* piv,
                        int* rank, const double* tol, double* work, int* info);

  // BLAS routine for Backward/Forward substitution
  void _FORTRAN(dtrsv)(const char * uplo, const char * trans, const char * diag, const int * n,
                       const double * a, const int * lda, double * x, const int * incx);
}

namespace Pita {

/* Constructor */

PivotedCholeskySolver::PivotedCholeskySolver(double tolerance) :
  RankDeficientSolver(tolerance), 
  choleskyFactor_()
{}

/* Mutators */

void
PivotedCholeskySolver::matrixIs(const SymFullMatrix & matrix) {
  int newSize = matrix.dim();
  choleskyFactor_.setSize(newSize);
  
  // Fill in lower triangular part of choleskyFactor
  const double * sourceData = const_cast<SymFullMatrix &>(matrix).data();
  for (int i = 0; i < newSize; ++i) {
    const double * sourceDataEnd = (sourceData + i) + 1;
    std::copy(sourceData, sourceDataEnd, choleskyFactor_[i]);
    sourceData = sourceDataEnd;
  }

  setMatrixSize(newSize);
  setFactorRank(0);
  setStatus(NON_FACTORIZED);
}

void
PivotedCholeskySolver::matrixIs(const FullSquareMatrix & matrix) {
  choleskyFactor_.copy(matrix);
  
  setMatrixSize(choleskyFactor_.dim());
  setFactorRank(0);
  setStatus(NON_FACTORIZED);
}

void
PivotedCholeskySolver::statusIs(Status s) {
  if (s == status())
    return;

  if (s == FACTORIZED) {
    // Perform factorization
    getFactorPermutation().sizeIs(matrixSize());

    if (matrixSize() > 0) {
      const char uplo = 'U';   // Lower triangular in C indexing == upper triangular in Fortran indexing 

      int rank;
      int info;
      SimpleBuffer<double> workspace(2 * matrixSize());

      _FORTRAN(dpstrf)(&uplo, &getMatrixSize(), choleskyFactor_.data(), &getMatrixSize(),
                       getFactorPermutation().array(), &rank, &getTolerance(), workspace.array(), &info);

      setFactorRank(rank);
    }
  }

  setStatus(s);
}

/* Read Accessor */

const Vector &
PivotedCholeskySolver::solution(Vector & rhs) const {
  if (rhs.size() != matrixSize()) {
    throw Fwk::RangeException("in PivotedCholeskySolver::solution - Size mismatch"); 
  }

  if (status() != FACTORIZED) {
    throw Fwk::RangeException("in PivotedCholeskySolver::solution - Non-factorized matrix");
  }

  // Pointer to the end of factorPermutation
  const int * fp_ptr_end = getFactorPermutation().array() + this->factorRank();
  const int * fp_ptr;

  // Permute rhs
  SimpleBuffer<double> vec(this->factorRank());
  fp_ptr = getFactorPermutation().array();
  for (double * vec_ptr = vec.array(); fp_ptr != fp_ptr_end; ++vec_ptr, ++fp_ptr) {
    *vec_ptr = rhs[*fp_ptr - 1]; // Offset to get C indexing from Fortran indexing
  }

  // Setup solve routine
  const char uplo = 'U';             // Lower triangular in C indexing == upper triangular in Fortran indexng
  const char forwardSubTrans = 'T';  // To solve Rt^{-1}
  const char backwardSubTrans = 'N'; // To solve R^{-1}
  const char diag = 'N';             // CholeskyFactor has non-unit diagonal
  const int incx = 1;                // Rhs vector elements are contiguous

  // Forward substitution
  _FORTRAN(dtrsv)(&uplo, &forwardSubTrans, &diag, &getFactorRank(), this->choleskyFactor().data(),
                  &getMatrixSize(), vec.array(), &incx);
  // Backward substitution
  _FORTRAN(dtrsv)(&uplo, &backwardSubTrans, &diag, &getFactorRank(), this->choleskyFactor().data(),
                  &getMatrixSize(), vec.array(), &incx);

  // Replace rhs by solution
  rhs.zero();
  fp_ptr = getFactorPermutation().array();
  for (const double * vec_ptr = vec.array(); fp_ptr != fp_ptr_end; ++vec_ptr, ++fp_ptr) {
    rhs[*fp_ptr - 1] = *vec_ptr; // Offset to get C indexing from Fortran indexing
  }

  return rhs;
}

} // end namespace Pita
