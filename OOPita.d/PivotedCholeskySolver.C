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
  
  // Fill in lower triangular part of choleskyFactor
  choleskyFactor_.setSize(newSize);
  const double * sourceData = const_cast<SymFullMatrix &>(matrix).data();
  for (int i = 0; i < newSize; ++i) {
    const double * sourceDataEnd = (sourceData + i) + 1;
    std::copy(sourceData, sourceDataEnd, choleskyFactor_[i]);
    sourceData = sourceDataEnd;
  }

  performFactorization();
}

void
PivotedCholeskySolver::transposedMatrixIs(const FullSquareMatrix & matrix) {
  choleskyFactor_.copy(matrix);
  performFactorization();
}

void
PivotedCholeskySolver::performFactorization() {
  setMatrixSize(choleskyFactor_.dim());
  
  setVectorSize(matrixSize());
  setOrdering(PERMUTED);
  getFactorPermutation().sizeIs(vectorSize());

  setFactorRank(0);
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


void
PivotedCholeskySolver::orderingIs(Ordering o) {
  if (ordering() == 0)
    return;

  if (o == COMPACT) {
    setVectorSize(factorRank());
    getFactorPermutation().sizeIs(vectorSize());

    for (int i = 0; i < factorRank(); ++i) {
      getFactorPermutation()[i] = i + 1; // Fortran numbering
    }
  }

  setOrdering(o);
}

/* Read Accessor */

const Vector &
PivotedCholeskySolver::solution(Vector & rhs) const {
  if (rhs.size() != vectorSize()) {
    throw Fwk::RangeException("in PivotedCholeskySolver::solution - Size mismatch"); 
  }

  // Permute rhs
  SimpleBuffer<double> perm_vec(this->factorRank());

  double * rhs_data;
  if (ordering() == PERMUTED) {
    // Pointers at the beginning / end of permutation array
    const int * fp_ptr = getFactorPermutation().array();
    const int * fp_ptr_end = getFactorPermutation().array() + this->factorRank();

    for (double * perm_vec_ptr = perm_vec.array(); fp_ptr != fp_ptr_end; ++perm_vec_ptr, ++fp_ptr) {
      *perm_vec_ptr = rhs[*fp_ptr - 1]; // Offset to get C indexing from Fortran indexing
    }
    rhs_data = perm_vec.array();
  } else {
    rhs_data = rhs.data();
  }

  // Setup solve routine
  const char uplo = 'U';             // Lower triangular in C indexing == upper triangular in Fortran indexng
  const char forwardSubTrans = 'T';  // To solve Rt^{-1}
  const char backwardSubTrans = 'N'; // To solve R^{-1}
  const char diag = 'N';             // CholeskyFactor has non-unit diagonal
  const int incx = 1;                // Rhs vector elements are contiguous

  // Forward substitution
  _FORTRAN(dtrsv)(&uplo, &forwardSubTrans, &diag, &getFactorRank(), this->choleskyFactor().data(),
                  &getMatrixSize(), rhs_data, &incx);
  // Backward substitution
  _FORTRAN(dtrsv)(&uplo, &backwardSubTrans, &diag, &getFactorRank(), this->choleskyFactor().data(),
                  &getMatrixSize(), rhs_data, &incx);

  if (ordering() == PERMUTED) {
    // Replace rhs by solution
    rhs.zero();
    
    // Pointers at the beginning / end of permutation array
    const int * fp_ptr = getFactorPermutation().array();
    const int * fp_ptr_end = getFactorPermutation().array() + this->factorRank();
    for (const double * perm_vec_ptr = perm_vec.array(); fp_ptr != fp_ptr_end; ++perm_vec_ptr, ++fp_ptr) {
      rhs[*fp_ptr - 1] = *perm_vec_ptr; // Offset to get C indexing from Fortran indexing
    }
  }

  return rhs;
}

} // end namespace Pita
