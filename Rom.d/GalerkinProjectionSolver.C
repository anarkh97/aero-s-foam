#include "GalerkinProjectionSolver.h"

#include <Utils.d/linkfc.h>

#include <complex>
#include <stdexcept>

extern "C" {
  // Cholesky factorization
  void _FORTRAN(dpotrf)(const char *uplo, const int *n,
                        double *a, const int *lda, int *info);

  // Solve with existing Cholesky factorization
  void _FORTRAN(dpotrs)(const char *uplo, const int *n, const int *nhrs,
                        const double *a, const int *lda,
                        double *b, const int *ldb, int *info);
}

template <>
void
GenGalerkinProjectionSolver<double>::factorReducedMatrix() {
  int info;
  _FORTRAN(dpotrf)("L", &basisSize_, reducedMatrix_.data(), &basisSize_, &info);

  assert(info == 0);
}

template <>
void
GenGalerkinProjectionSolver<double>::solveReducedRhs() {
  const int INT_ONE = 1;
  
  int info;
  _FORTRAN(dpotrs)("L", &basisSize_, &INT_ONE,
                   reducedMatrix_.data(), &basisSize_,
                   reducedRhs_.data(), &basisSize_,
                   &info);
  
  assert(info == 0);
}

template <>
void
GenGalerkinProjectionSolver<std::complex<double> >::factorReducedMatrix() {
  throw std::logic_error("Not implemented");
}

template <>
void
GenGalerkinProjectionSolver<std::complex<double> >::solveReducedRhs() {
  throw std::logic_error("Not implemented");
}
