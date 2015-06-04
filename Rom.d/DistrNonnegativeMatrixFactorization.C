#if defined(USE_SCALAPACK) && defined(USE_EIGEN3)
#include "DistrNonnegativeMatrixFactorization.h"

#include <Comm.d/Communicator.h>
#include <Rom.d/SparseSolvers.d/ScalaLH.d/Plh.h>
#include <Utils.d/linkfc.h>

#include <mpi.h>

extern "C" {
  // Context & cpu topology management
  int Csys2blacs_handle(MPI_Comm comm);

  void Cfree_blacs_system_handle(int handle);
  void Cblacs_gridinit(int *ictxt, char *order, int nprow, int npcol);

  void Cblacs_gridinfo(int ictxt, int *nprow, int *npcol, int *myrow, int *mycol);
  void Cblacs_gridexit(int ictxt);

  // Index mapping 
  int _FORTRAN(indxg2l)(const int *indxglob, const int *nb, const int *iproc_dummy, const int *isrcproc_dummy, const int *nprocs);
  int _FORTRAN(indxg2p)(const int *indxglob, const int *nb, const int *iproc_dummy, const int *isrcproc, const int *nprocs);
}

namespace Rom {

DistrNonnegativeMatrixFactorization
::DistrNonnegativeMatrixFactorization(Communicator * comm, int rowCount, int colCount, int localRows, int basisDimension,
                                      int blockSize, int maxIter, double tol, int method) : 
  communicator_(comm),
  rowCount_(rowCount),
  colCount_(colCount),
  localRows_(localRows),
  basisDimension_(basisDimension),
  blockSize_(blockSize),
  maxIter_(maxIter),
  tol_(tol),
  method_(method),
  matrixBuffer_(localRows,colCount),
  basisBuffer_(localRows,basisDimension)
{
}

void
DistrNonnegativeMatrixFactorization::solve()
{
  // 1. Construct matrix A
  int blacsHandle = Csys2blacs_handle(*communicator_->getCommunicator());
  int context = blacsHandle;
  char order[] = "R";
  int rowCpus = communicator_->numCPUs();
  int colCpus = 1;
  Cblacs_gridinit(&context, order, rowCpus, colCpus);
  SCDoubleMatrix A(context, rowCount_, colCount_, blockSize_, blockSize_);

  // 2. Copy matrixBuffer into A
  int myrow, mycol;
  int dummy, zero=0;
  Cblacs_gridinfo(context, &rowCpus, &colCpus, &myrow, &mycol);
  double *row = new double[std::max(colCount_,basisDimension_)];
  for(int j=1; j<=rowCount_; j++) {
    int p = _FORTRAN(indxg2p)(&j, &blockSize_, &dummy, &zero, &rowCpus);
    if(myrow == p) {
      int lj = _FORTRAN(indxg2l)(&j, &blockSize_, NULL, NULL, &rowCpus) - 1;
      for(int k=0; k<colCount_; ++k) row[k] = -matrixBuffer_(lj,k);
      A.setMatrixRow(j, row);
    }
  }

  // Construct matrices to store factors W and H
  SCDoubleMatrix W(context, rowCount_, basisDimension_, blockSize_, blockSize_);
  SCDoubleMatrix Htranspose(context, colCount_, basisDimension_, blockSize_, blockSize_);

  // Initialize W with random positive entries between 0 and 1
  W.initRandom(1,0.,1.);

  switch(method_) {
    default: case 1 : {
      // Alternating non-negative least squares iteration loop
      for(int i=0; i<maxIter_; ++i) {

        // make a copy of W to use for checking the stopping criteria
        SCDoubleMatrix W_copy(W);

        // solve: min ||WH-A|| s.t. H >= 0
        solveNNLS_MRHS(W, A, Htranspose, 0);

        // solve: min ||H^TW^T-A^T|| s.t. W >= 0
        solveNNLS_MRHS(Htranspose, A, W, 1);

        // compute residual and check stopping criteria
        SCDoubleMatrix Err(A); W.multiply(Htranspose, Err, 'N', 'T', -1.0, 1.0); // Err = A-W*H;
        double res = Err.froNorm()/A.froNorm();
        W.add(W_copy, 'N', rowCount_, basisDimension_, 1.0, -1.0); // W_copy = W-W_copy
        double inc = W_copy.froNorm();
        if(communicator_->myID() == 0)
          std::cout << "iteration = " << i+1 << ", rel. residual = " << res << ", solution incr. = " << inc << std::endl;
        if(inc < tol_) break;
      }
    } break;

    case 2 : {
      if(communicator_->myID() == 0)
        std::cout << "ERROR: Greedy method is not implemented in DistrNonnegativeMatrixFactorization.C\n";
      exit(-1);
    } break;
  }

  // copy W into basisBuffer_
  for(int j=1; j<=rowCount_; j++) {
    int p = _FORTRAN(indxg2p)(&j, &blockSize_, &dummy, &zero, &rowCpus);
    if(myrow == p) {
      int lj = _FORTRAN(indxg2l)(&j, &blockSize_, NULL, NULL, &rowCpus) - 1;
      W.getMatrixRow(j, row, 'U');
      for(int k=0; k<basisDimension_; ++k) basisBuffer_(lj,k) = row[k];
    }
  }

  delete [] row;

  Cblacs_gridexit(context);
  Cfree_blacs_system_handle(blacsHandle);
}

void
DistrNonnegativeMatrixFactorization::solveNNLS_MRHS(SCDoubleMatrix &A, SCDoubleMatrix &B, SCDoubleMatrix &X, int flag)
{
  // if flag = 0 then solve: min ||AX^T-B|| s.t. X >= 0
  // if flag = 1 then solve: min ||AX^T-B^T|| s.t. X >= 0

  Plh solver(A.getNumberOfRows(), A.getNumberOfCols());
  solver.init();
  solver.setMatrix(A);
  solver.setVerbose(0);

  double *rhs = new double[A.getNumberOfRows()];
  double *sol = new double[A.getNumberOfCols()];
  int nrhs = X.getNumberOfRows();

  for(int i=1; i<=nrhs; ++i) {
    if(flag==0) B.getMatrixColumn(i,rhs,'A');
    else        B.getMatrixRow(i,rhs,'A');
    solver.setRHS(rhs);
    solver.setRtol(1e-16);
    solver.solve();
    solver.getSolution(sol);
    X.setMatrixRow(i,sol);
  }

  delete [] rhs;
  delete [] sol;
}

} // end namespace Rom
#endif
