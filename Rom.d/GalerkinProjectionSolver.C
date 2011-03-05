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

namespace Rom {

template <>
void
GenGalerkinProjectionSolver<double>::performFactor() {
  const int basisDim = basisSize();

  int info;
  _FORTRAN(dpotrf)("L", &basisDim, reducedMatrix_.data(), &basisDim, &info);

  assert(info == 0);
}

template <>
void
GenGalerkinProjectionSolver<double>::performSolve() {
  const int basisDim = basisSize();
  const int INT_ONE = 1;
  
  int info;
  _FORTRAN(dpotrs)("L", &basisDim, &INT_ONE,
                   reducedMatrix_.data(), &basisDim,
                   getReducedSolution().data(), &basisDim,
                   &info);
  
  assert(info == 0);
}

template <>
void
GenGalerkinProjectionSolver<std::complex<double> >::performFactor() {
  throw std::logic_error("Not implemented");
}

template <>
void
GenGalerkinProjectionSolver<std::complex<double> >::performSolve() {
  throw std::logic_error("Not implemented");
}

} /* end namespace Rom */
