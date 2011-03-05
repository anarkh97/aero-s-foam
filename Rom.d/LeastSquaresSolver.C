#include "LeastSquaresSolver.h"

#include <Utils.d/linkfc.h>

#include <complex>

#include <cassert>

extern "C" {
  // Perform the QR factorization
  void _FORTRAN(dgeqrf)(const int *m, const int *n, double *a, const int *lda,
                        double *tau, double *work, const int *lwork, int *info);

  // Solves Q^{-1} = Q^T
  void _FORTRAN(dormqr)(const char *side, const char *trans,
                        const int *m, const int *n, const int *k,
                        const double *a, const int *lda, const double *tau,
                        double *c, const int *ldc,
                        double *work, const int *lwork, int *info);

  // Solves R^{-1}
  void _FORTRAN(dtrtrs)(const char *uplo, const char *trans, const char* diag,
                        const int *n, const int *nrhs,
                        const double *a, const int *lda,
                        double *b, const int *ldb, int *info);
  
  // Multiplies by R
  void _FORTRAN(dtrmm)(const char *side, const char *uplo, const char *transa, const char* diag,
                       const int *m, const int *n,
                       const double *alpha, const double *a, const int *lda,
                       double *b, const int *ldb);
}

namespace Rom {

// Implementation specializations

template <>
void
GenLeastSquaresSolver<double>::factor() {
  int info;

  const int lworkQuery = -1;
  double lworkAns;
  _FORTRAN(dgeqrf)(&equationCount_, &unknownCount_, matrixBuffer_.array(), &matrixLeadDim_,
                   tauBuffer_.array(), &lworkAns, &lworkQuery, &info);
  assert(info == 0);

  const int lwork = static_cast<int>(lworkAns);
  ScalarBuffer workBuffer(lwork);
  _FORTRAN(dgeqrf)(&equationCount_, &unknownCount_, matrixBuffer_.array(), &matrixLeadDim_,
                   tauBuffer_.array(), workBuffer.array(), &lwork, &info);
  assert(info == 0);
}

template <>
void
GenLeastSquaresSolver<double>::project() {
  const char left = 'L';
  const char trans = 'T';

  int info;
  
  const int lworkQuery = -1;
  double lworkAns;
  _FORTRAN(dormqr)(&left, &trans,
                   &equationCount_, &rhsCount_, &unknownCount_,
                   matrixBuffer_.array(), &matrixLeadDim_, tauBuffer_.array(), 
                   rhsBuffer_.array(), &rhsLeadDim_,
                   &lworkAns, &lworkQuery, &info);
  assert(info == 0);

  const int lwork = static_cast<int>(lworkAns);
  ScalarBuffer workBuffer(lwork);
  _FORTRAN(dormqr)(&left, &trans,
                   &equationCount_, &rhsCount_, &unknownCount_,
                   matrixBuffer_.array(), &matrixLeadDim_, tauBuffer_.array(), 
                   rhsBuffer_.array(), &rhsLeadDim_,
                   workBuffer.array(), &lwork, &info);
  assert(info == 0);
}
 
template <>
void
GenLeastSquaresSolver<double>::solve() {
  const char upper = 'U';
  const char no = 'N'; // not-transposed, non-unit

  int info;
  
  _FORTRAN(dtrtrs)(&upper, &no, &no,
                   &unknownCount_, &rhsCount_,
                   matrixBuffer_.array(), &matrixLeadDim_,
                   rhsBuffer_.array(), &rhsLeadDim_, &info);
  assert(info == 0);
}

template <>
void
GenLeastSquaresSolver<double>::unsolve() {
  const char left = 'L';
  const char upper = 'U';
  const char no = 'N'; // not-transposed, non-unit
  const double one = 1.0;

  _FORTRAN(dtrmm)(&left, &upper, &no, &no,
                  &unknownCount_, &rhsCount_,
                  &one, matrixBuffer_.array(), &matrixLeadDim_,
                  rhsBuffer_.array(), &rhsLeadDim_);
}

// Interface specializations

using LeastSquares::READY;
using LeastSquares::FACTORED;
using LeastSquares::PROJECTED;
using LeastSquares::SOLVED;

template <>
void
GenLeastSquaresSolver<double>::statusIs(Status s) {
  if (s > status()) {
    while (status() != s) {
      switch (status()) {
        case READY:
          factor(); break;
        case FACTORED:
          project(); break;
        case PROJECTED:
          solve(); break;
        default:
          throw std::logic_error("Internal error");
      }

      status_ = Status(status_ + 1);
    }
  } else {
    if (status() == SOLVED && s == PROJECTED) {
      unsolve();
    }

    status_ = s;
  }
}

template <>
void
GenLeastSquaresSolver<std::complex<double> >::statusIs(Status) {
  std::logic_error("Not implemented");
}

} /* end namespace Rom */
