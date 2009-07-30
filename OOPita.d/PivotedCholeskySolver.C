#include "PivotedCholeskySolver.h"

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

PivotedCholeskySolver::PivotedCholeskySolver() :
  status_(NON_FACTORIZED),
  matrixSize_(0),
  factorRank_(0),
  choleskyFactor_(),
  factorPermutation_()
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

  // Remark: factorRank_, factorPermutation_ not up-to-date
  matrixSize_ = newSize;
  status_ = NON_FACTORIZED;
}

void
PivotedCholeskySolver::matrixIs(const FullSquareMatrix & matrix) {
  choleskyFactor_.copy(matrix);
  
  // Remark: factorRank_, factorPermutation_ not up-to-date
  matrixSize_ = choleskyFactor_.dim();
  status_ = NON_FACTORIZED;
}

void PivotedCholeskySolver::statusIs(Status s) {
  if (s == status())
    return;

  if (s == FACTORIZED) {
    // Perform factorization
    factorPermutation_.sizeIs(matrixSize());

    if (matrixSize() > 0) {
      const char uplo = 'U';   // Lower triangular in C indexing == upper triangular in Fortran indexing 
      const double tol = -1.0; // Use default tolerance threshold

      int info;
      SimpleBuffer<double> workspace(2 * matrixSize_);

      _FORTRAN(dpstrf)(&uplo, &matrixSize_, choleskyFactor_.data(), &matrixSize_,
                       factorPermutation_.array(), &factorRank_, &tol, workspace.array(), &info);
    }
  }

  status_ = s;
}

/* Read Accessor */

Vector &
PivotedCholeskySolver::solution(Vector & rhs) const {
  if (rhs.size() != this->matrixSize()) {
    throw Fwk::RangeException("in PivotedCholeskySolver::solution - Size mismatch"); 
  }

  if (this->status() != FACTORIZED) {
    throw Fwk::RangeException("in PivotedCholeskySolver::solution - Non-factorized matrix");
  }

  // Pointer to the end of factorPermutation
  const int * fp_ptr_end = factorPermutation_.array() + this->factorRank();
  const int * fp_ptr;

  // Permute rhs
  SimpleBuffer<double> vec(this->factorRank());
  fp_ptr = factorPermutation_.array();
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
  _FORTRAN(dtrsv)(&uplo, &forwardSubTrans, &diag, &factorRank_, this->choleskyFactor().data(),
                  &matrixSize_, vec.array(), &incx);
  // Backward substitution
  _FORTRAN(dtrsv)(&uplo, &backwardSubTrans, &diag, &factorRank_, this->choleskyFactor().data(),
                  &matrixSize_, vec.array(), &incx);

  // Replace rhs by solution
  rhs.zero();
  fp_ptr = factorPermutation_.array();
  for (const double * vec_ptr = vec.array(); fp_ptr != fp_ptr_end; ++vec_ptr, ++fp_ptr) {
    rhs[*fp_ptr - 1] = *vec_ptr; // Offset to get C indexing from Fortran indexing
  }

  return rhs;
}

} // end namespace Pita
