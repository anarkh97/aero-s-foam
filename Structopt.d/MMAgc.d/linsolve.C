#ifdef STRUCTOPT

#include <stdlib.h>
#include <stdio.h>
#include <Utils.d/dbg_alloca.h>
#include <Utils.d/linkfc.h>

#include "mma.h"

extern "C" {

   void  _FORTRAN(dgesv) (int&, int&, double*, int&, int*, 
                          double*, int&, int&);

} ;

void MMAgc::linsolve(double **A, double* x, double *b, int n) {

  // interface routine to LAPACK function: dgesv
  // 
  // solves  A x = b
  //
  // where A : (n x n), symmetric
  //       x : (n x 1)
  //       b : (n x 1)
  // 
  // watch: A will be factorized
  // 

  int nrhs = 1;
  int lda  = n;
  int ldb  = n;
  int info;  

  int* ipiv = static_cast<int*>(dbg_alloca(n*sizeof(int)));
  
  int i;
  for (i=0;i<n;i++) x[i] = b[i];

  _FORTRAN(dgesv) (n, nrhs, A[0], lda, ipiv, x, ldb, info);

  if (info > 0) {
    fprintf(stderr,"Error in linsolve: DGESV returned %d\n",info);
    exit(-1);
  }
}
  
#endif
