#include "CholeskyUtils.h"

#include <Math.d/FullSquareMatrix.h>

#include <Utils.d/linkfc.h>

#include <complex>
#include <stdexcept>

extern "C" {
  // Cholesky factorization
  void _FORTRAN(dpotrf)(const char *uplo, const int *n,
                        double *a, const int *lda, int *info);

  // Inverse triangular matrix
  void _FORTRAN(dtrtri)(const char *uplo, const char *diag, const int *n,
                        double *a, const int *lda, int *info);
}

namespace Rom {

template <>
const GenFullSquareMatrix<double> &
cholesky_factor_upper(GenFullSquareMatrix<double> &m) {
  const int basisDim = m.dim();

  int info;
  _FORTRAN(dpotrf)("L", &basisDim, m.data(), &basisDim, &info);

  if (info) {
    throw std::runtime_error("Error in dpotrf (Cholesky factorization)");
  };

  return m;
}

template <>
const GenFullSquareMatrix<double> &
cholesky_factor_lower(GenFullSquareMatrix<double> &m) {
  const int basisDim = m.dim();

  int info;
  _FORTRAN(dpotrf)("U", &basisDim, m.data(), &basisDim, &info);

  if (info) {
    throw std::runtime_error("Error in dpotrf (Cholesky factorization)");
  };

  return m;
}

template <>
const GenFullSquareMatrix<double> &
inverse_triangular_upper(GenFullSquareMatrix<double> &m) {
  const int basisDim = m.dim();

  int info;
  _FORTRAN(dtrtri)("L", "N", &basisDim, m.data(), &basisDim, &info);

  if (info) {
    throw std::runtime_error("Error in dtrtri (Triangular matrix inversion)");
  };

  return m;
}


template <>
const GenFullSquareMatrix<double> &
inverse_triangular_lower(GenFullSquareMatrix<double> &m) {
  const int basisDim = m.dim();

  int info;
  _FORTRAN(dtrtri)("U", "N", &basisDim, m.data(), &basisDim, &info);

  if (info) {
    throw std::runtime_error("Error in dtrtri (Triangular matrix inversion)");
  };

  return m;
}

//-------------------------------------------------------------------

template <>
const GenFullSquareMatrix<std::complex<double> > &
cholesky_factor_upper(GenFullSquareMatrix<std::complex<double> > &) {
  throw std::logic_error("Not implemented");
}

template <>
const GenFullSquareMatrix<std::complex<double> > &
cholesky_factor_lower(GenFullSquareMatrix<std::complex<double> > &) {
  throw std::logic_error("Not implemented");
}

template <>
const GenFullSquareMatrix<std::complex<double> > &
inverse_triangular_upper(GenFullSquareMatrix<std::complex<double> > &) {
  throw std::logic_error("Not implemented");
}

template <>
const GenFullSquareMatrix<std::complex<double> > &
inverse_triangular_lower(GenFullSquareMatrix<std::complex<double> > &) {
  throw std::logic_error("Not implemented");
}

} /* end namespace Rom */
