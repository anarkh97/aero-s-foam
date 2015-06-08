#if defined(USE_SCALAPACK) && defined(USE_EIGEN3)
#include "DistrNonnegativeMatrixFactorization.h"

#include <Comm.d/Communicator.h>
#include <Rom.d/SparseSolvers.d/ScalaLH.d/Plh.h>
#include <Utils.d/linkfc.h>
#include <Timers.d/GetTime.h>

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

extern int verboseFlag;

namespace Rom {

double t1=0,t2=0,t3=0,t4=0,t5=0,t6=0;

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
  SCDoubleMatrix A(context, rowCount_, colCount_, blockSize_, blockSize_, *communicator_->getCommunicator());

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
  SCDoubleMatrix W(context, rowCount_, basisDimension_, blockSize_, blockSize_, *communicator_->getCommunicator());
  SCDoubleMatrix Htranspose(context, colCount_, basisDimension_, blockSize_, blockSize_, *communicator_->getCommunicator());

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

  for(int i=0; i<communicator_->numCPUs(); ++i) { if(i==communicator_->myID()) {
  std::cerr << "t1 = " << t1/1000 << ", t2 = " << t2/1000 << ", t3 = " << t3/1000 << ", t4 = " << t4/1000 
            << ", t5 = " << t5/1000 << ", t6 = " << t6/1000 << std::endl;
  } communicator_->sync(); }
}

void
DistrNonnegativeMatrixFactorization::solveNNLS_MRHS(SCDoubleMatrix &A, SCDoubleMatrix &B, SCDoubleMatrix &X, int flag)
{
  // if flag = 0 then solve: min ||AX^T-B|| s.t. X >= 0
  // if flag = 1 then solve: min ||AX^T-B^T|| s.t. X >= 0
  if(verboseFlag && communicator_->myID() == 0) {
    std::cerr << "\nm = " << A.getNumberOfRows() 
              << ", n = " << A.getNumberOfCols()
              << ", k = " << X.getNumberOfRows() << std::endl;
  }

  // number of sub-matrices in column-wise partition of B if flag == 0, or row-wise partition of B if flag == 1
  const int nsub = communicator_->numCPUs(); // XXX

  // create new communicator and context if necessary
  MPI_Comm comm;
  int color, blacsHandle, context, rowCpus, colCpus;
  if(nsub > 1) {
    color = communicator_->myID()%nsub;
    MPI_Comm_split(*communicator_->getCommunicator(), color+1, 0, &comm);
    blacsHandle = Csys2blacs_handle(comm);
    context = blacsHandle;
    char order[] = "R";
    MPI_Comm_size(comm, &rowCpus);
    colCpus = 1;
    Cblacs_gridinit(&context, order, rowCpus, colCpus);
  }
  else {
    color = 0;
    comm = *communicator_->getCommunicator();
    context = A.getContext();
    rowCpus = A.getNumberOfProcsRow();
    colCpus = A.getNumberOfProcsCol();
  }

  Plh solver(A.getNumberOfRows(), A.getNumberOfCols());
  solver.setContext(context, rowCpus, colCpus, comm);
  solver.setQProcGrid(rowCpus, colCpus);
  solver.setABlockSize(blockSize_, blockSize_);
  solver.setQBlockSize(blockSize_, blockSize_);
  //solver.setColumnScaling();
  solver.init();
  t1 -= getTime();
  for(int k = 0; k < nsub; ++k) {
    if(k != color) A.copyRedist(A.getNumberOfRows(), A.getNumberOfCols(), 1, 1, A.getContext());
    else solver.setMatrix(A);
  }
  t1 += getTime();
  solver.setVerbose(0);
  solver.setMaxIterRatio(3);

  SCDoubleMatrix &b = solver.getRhsVector(); 
  SCDoubleMatrix &x = solver.getSolutionVector();

  // copy/redistribute from B to subB, if necessary
  SCDoubleMatrix *subB, *subX;
  int nrhs;
  if(nsub == 0) {
    nrhs = X.getNumberOfRows();
    subX = &X;
    subB = &B;
  }
  else {
    nrhs = X.getNumberOfRows()/nsub + ((color < X.getNumberOfRows()%nsub) ? 1 : 0);
    subX = new SCDoubleMatrix(context, nrhs, X.getNumberOfCols(), blockSize_, blockSize_, comm);
    t5 -= getTime();
    if(flag == 0) {
      subB = new SCDoubleMatrix(context, B.getNumberOfRows(), nrhs, blockSize_, blockSize_, comm);
      for(int k=0,ja=1,n; k<nsub; ++k, ja+=n) {
        n = X.getNumberOfRows()/nsub + ((k < X.getNumberOfRows()%nsub) ? 1 : 0);
        if(k != color) B.copyRedist(A.getNumberOfRows(), n, 1, ja, A.getContext());              // send only
        else           B.copyRedist(A.getNumberOfRows(), n, 1, ja, *subB, 1, 1, A.getContext()); // send and receive
      }
    }
    else {
      subB = new SCDoubleMatrix(context, nrhs, B.getNumberOfCols(), blockSize_, blockSize_, comm);
      for(int k=0,ia=1,m; k<nsub; ++k, ia+=m) {
        m = X.getNumberOfRows()/nsub + ((k < X.getNumberOfRows()%nsub) ? 1 : 0);
        if(k != color) B.copyRedist(m, A.getNumberOfRows(), ia, 1, A.getContext());              // send only
        else           B.copyRedist(m, A.getNumberOfRows(), ia, 1, *subB, 1, 1, A.getContext()); // send and receive
      }
    }
    t5 += getTime();
  }

  int iter = 0;
  for(int i=1; i<=nrhs; i++) {
    
    t2 -= getTime();
    // copy ith column/row of matrix subB into vector b
    if(flag==0) subB->add(b, 'N', A.getNumberOfRows(), 1, 1.0, 0.0, 1, i, 1, 1);
    else        subB->add(b, 'T', A.getNumberOfRows(), 1, 1.0, 0.0, i, 1, 1, 1);
    t2 += getTime();

    // solve: min ||Ax-b|| s.t. x >= 0
    t3 -= getTime();
    solver.setRtol(1e-16);
    solver.solve();
    iter += solver.getIter();
    t3 += getTime();

    // copy x to ith row of subX
    t4 -= getTime();
    x.add(*subX, 'N', 1, A.getNumberOfCols(), 1.0, 0.0, 1, 1, i, 1);
    t4 += getTime();
  }

  if(verboseFlag && communicator_->myID() == 0) {
    std::cerr << "Total number of iterations = " << iter << std::endl;
  }
  //solver.printTimes(true);

  // copy/redistribute from subX to X, if necessary
  if(nsub > 1) {
    t6 -= getTime();
    for(int k=0,ib=1,m; k<nsub; ++k,ib+=m) {
      m = X.getNumberOfRows()/nsub + ((k < X.getNumberOfRows()%nsub) ? 1 : 0);
      if(k != color) SCDoubleMatrix::copyRedist(m, A.getNumberOfCols(), X, ib, 1, A.getContext()); // receive only
      else           subX->copyRedist(m, A.getNumberOfCols(), 1, 1, X, ib, 1, A.getContext());     // send and receive
    }
    t6 += getTime();
    delete subB;
    delete subX;
    Cblacs_gridexit(context);
    Cfree_blacs_system_handle(blacsHandle);
  }
}

} // end namespace Rom
#endif
