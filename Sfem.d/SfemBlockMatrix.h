#ifndef SFEM_BLOCK_MATRIX_H_
#define SFEM_BLOCK_MATRIX_H_

#include<Utils.d/MyComplex.h>
#include<Math.d/matrix.h>

template<class Scalar>
class SfemBlockMatrix {
  int n, L, P, ndim;   // L & P comes from class Sfem, n comes from domain
  Scalar*** kiuj;
  GenSparseMatrix<Scalar>** allK;
  double ***inpsi;
  Connectivity *blockToBlock;
  Scalar *diags;
  int *firstdof;
  GenFullM<Scalar> ** diagblocks;
  GenSolver<Scalar>* meansolver;
  Scalar* scalarfactors;
  BlockInfo inf; 
 public:
  SfemBlockMatrix(int L, int n, int P, int _ndim, int output_order);
  void setKi(GenSparseMatrix<Scalar>* K, int i);
  void matvec_sfem_block(GenVector<Scalar> &u);
  void mult(GenVector<Scalar> &u, GenVector<Scalar> &ku);
  BlockInfo & dim(); 
  Scalar diag(int i);
  int* getFirstDof();
  int numNodes();
  GenFullM<Scalar>* getDiagMatrix(int i);
  void setMeanSolver(GenSolver<Scalar> *prc);
  GenSolver<Scalar>* getMeanSolver();
  Scalar* getBlockScalarMultipliers();
  int getBlockSize() {return n;}
  int neqs() {return n*P;}
};

#ifdef _TEMPLATE_FIX_
#include <Sfem.d/SfemBlockMatrix.C>
#endif

#endif
