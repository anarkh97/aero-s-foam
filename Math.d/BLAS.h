//
// Created by Michel Lesoinne on 11/14/17.
//

#ifndef FEM_BLAS_H
#define FEM_BLAS_H
#include <complex>
#include <Utils.d/linkfc.h>
extern "C" {

  void _FORTRAN(dgemm)(const char &, const char &, const int &,const int &,
                       const int &, const double &, double *, const int &,
                       double *, const int &, const double &, double *, const int &);

  void _FORTRAN(zgemm)(const char &, const char &, const int &,const int &,
                       const int &, const std::complex<double> &, std::complex<double> *, const int &,
                       std::complex<double> *, const int &, const std::complex<double> &, std::complex<double> *,
                       const int &);

  void _FORTRAN(dgemv)(const char &, const int &,const int &,
                       const double &, const double *, const int &,
                       const double *, const int &, const double &, double *, const int &);

  void _FORTRAN(zgemv)(const char &, const int &,const int &,
                       const std::complex<double> &, std::complex<double> *, const int &,
                       std::complex<double> *, const int &, const std::complex<double> &, std::complex<double> *, const int &);

  void _FORTRAN(dgesv)(const int &, const int &, double *, const int &, int *, double *,
                       const int &, int &);

  void _FORTRAN(zgesv)(const int &, const int &, std::complex<double> *, const int &, int *, std::complex<double> *,
                       const int &, int &);

}
#endif //FEM_BLAS_H
