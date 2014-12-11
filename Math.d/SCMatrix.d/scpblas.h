#ifndef SCPBLAS_H_
#define SCPBLAS_H_

#ifdef NNLS_DEV
#include "linkfc.h"
#else
#include <Utils.d/linkfc.h>
#endif

extern "C" {
#ifdef USE_MKL
#include <mkl_pblas.h>
#include <mkl_blacs.h>
#else
#include <pblas.h>
#include <PBpblas.h>
#include <PBblacs.h>
#endif
}

// Need to replace this with an appropriate blacs header file.
// Adjusted for 0 based indexing
#define DLEN_ 9
#define DTYPE_ 0
#define CTXT_ 1
#define M_ 2
#define N_ 3
#define MB_ 4
#define NB_ 5
#define RSRC_ 6
#define CSRC_ 7
#define LLD_ 8


#ifdef __cplusplus
extern "C" {
    int  _FORTRAN(numroc)(int *n, int *nb, int *myrow, int *zero, int *nprow );
    void _FORTRAN(descinit)(int *descb, int *n, int *uno, int *mb, int *nb, int *mm, int *nn, int *context, int *lld, int *info );
    void _FORTRAN(pdelset)( double *b, int *i, int *uno, int *descb, double *alpha );
    void _FORTRAN(pdelget)( char *, char *, double *, double *, int *, int *, int *);
    void _FORTRAN(pielset)( int *b, int *i, int *uno, int *descb, int *alpha );
    void _FORTRAN(pielget)( char *, char *, int *, int *, int *, int *, int *);
    void _FORTRAN(pdgesv)(int *, int *, double *, int *, int *, int *, int *, double *, int *, int *, int *, int *);
    int  _FORTRAN(indxl2g)(int *, int *, int *, int *, int *);
    int  _FORTRAN(indxg2l)(int *, int *, int *, int *, int *);
    int  _FORTRAN(indxg2p)(int *, int *, int *, int *, int *);
    void _FORTRAN(pdlaprnt)( int *, int *, double *, int *, int *, int *, int *, int *, char**, int *, double *);
    void _FORTRAN(pdlapiv)(char *, char *, char *, int *, int *, double *, int *, int *, int *, int *, int *, int *, int *, int *);
    void _FORTRAN(fpdlapiv)(char *, char *, char *, int *, int *, double *, int *, int *, int *, int *, int *, int *, int *, int *);
    void _FORTRAN(fpslapiv)(char *, char *, char *, int *, int *, float *, int *, int *, int *, int *, int *, int *, int *, int *);
    void _FORTRAN(infog2l)(int *, int *, int*, int *, int *, int*, int *, int *, int*, int*, int*);
    void _FORTRAN(pdgeqrf)(int *, int *, double *, int *, int *, int *, double *, double *, int *, int *);
    void _FORTRAN(pdgels)(char *, int *, int *, int *, double *, int *, int *, int *,
                          double *, int *, int *, int *, double *, int *, int *);
    void _FORTRAN(pdormqr)(char *, char *, int *, int *, int *, double *, int *, int *, int *, double *,
                           double *, int *, int *, int *, double *, int *, int *);
    void _FORTRAN(pdtrtrs)(char *, char *, char *, int *, int *, double *, int *, int *, int *,
                           double *, int *, int *, int *, int *);
    void _FORTRAN(pdlacp2)(char *, int *, int *, double *, int *, int *, int *, double *, int *, int *, int *);
    int  _FORTRAN(ilcm)(int *, int *);
    int  _FORTRAN(pdgemr2d)(int *, int *, double *, int *, int *, int *, double *, int *, int *, int *, int *);
    double _FORTRAN(pdlange)(char *, int *, int *, double *, int *, int *, int *, double *);
}
#endif

#endif /* SCPBLAS_H_ */
