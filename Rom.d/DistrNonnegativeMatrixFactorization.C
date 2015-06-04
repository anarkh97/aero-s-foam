#include "DistrNonnegativeMatrixFactorization.h"

#include <Comm.d/Communicator.h>
#include <Rom.d/SparseSolvers.d/ScalaLH.d/Plh.h>
#include <Utils.d/linkfc.h>

#ifdef USE_SCALAPACK
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

#endif /* USE_SCALAPACK */

namespace Rom {

DistrNonnegativeMatrixFactorization
::DistrNonnegativeMatrixFactorization(Communicator * comm, int rowCount, int colCount, int localRows, int basisDimension, int blockSize, int maxIter, double tol) : 
  communicator_(comm),
  rowCount_(rowCount),
  colCount_(colCount),
  localRows_(localRows),
  basisDimension_(basisDimension),
  blockSize_(blockSize),
  maxIter_(maxIter),
  tol_(tol),
  matrixBuffer_(localRows,colCount),
  basisBuffer_(localRows,basisDimension)
{
}

DistrNonnegativeMatrixFactorization::~DistrNonnegativeMatrixFactorization()
{
}

void
DistrNonnegativeMatrixFactorization::solve(int)
{
#ifdef USE_SCALAPACK
  int blacsHandle = Csys2blacs_handle(*communicator_->getCommunicator());

  int context = blacsHandle;
  char order[] = "R";
  int rowCpus = communicator_->numCPUs();
  int colCpus = 1;
  Cblacs_gridinit(&context, order, rowCpus, colCpus);

  SCDoubleMatrix A(context, rowCount_, colCount_, blockSize_, blockSize_);

  int myrow, mycol;
  int dummy, zero=0;
  Cblacs_gridinfo(context, &rowCpus, &colCpus, &myrow, &mycol);
  double *row = new double[std::max(colCount_,basisDimension_)];
  for(int i=0; i<communicator_->numCPUs(); ++i) { if(i==communicator_->myID()) {
  std::cerr << "rowCount_ = " << rowCount_ << ", colCount_ = " << colCount_ << ", localRows_ = " << localRows_ << std::endl;
  for(int j=1; j<=rowCount_; j++) {
    int p = _FORTRAN(indxg2p)(&j, &blockSize_, &dummy, &zero, &rowCpus);
    if(myrow == p) {
      int lj = _FORTRAN(indxg2l)(&j, &blockSize_, NULL, NULL, &rowCpus) - 1;
      std::cout << "j = " << j << ", p = " << p << ", lj = " << lj << std::endl;
      for(int k=0; k<colCount_; ++k) row[k] = -matrixBuffer_(lj,k);
      A.setMatrixRow(j, row);
    }
  }
  } communicator_->sync(); }

  double Anorm = A.froNorm();
  for(int i=0; i<communicator_->numCPUs(); ++i) { if(i==communicator_->myID()) {
  std::cerr << "norm of A = " << Anorm << std::endl;
  } communicator_->sync(); }

  //basisBuffer_.setZero();

  SCDoubleMatrix W(context, rowCount_, basisDimension_, blockSize_, blockSize_);
  SCDoubleMatrix Htranspose(context, colCount_, basisDimension_, blockSize_, blockSize_);

  // initialize W with random positive entries
  W.initRandom(); // XXX W.setRandom(); W += Eigen::MatrixXd::Ones(m,k);

  SCDoubleMatrix W_copy(W);
  for(int i=0; i<maxIter_; ++i) {
    solveNNLS_MRHS(W, A, Htranspose, 0); // solve: min ||WH-A|| s.t. H >= 0
    solveNNLS_MRHS(Htranspose, A, W, 1); // solve: min ||H^TW^T-A^T|| s.t. W >= 0

    // compute residual and check convergence criteria
/* XXX
    SCDoubleMatrix Err(A); Err.multiply(W,H,1.0,-1.0); // Err = A-W*H;
    int index = findColumnWithLargestMagnitude(Err);
    res = Err.froNorm()/A.froNorm();
    double inc = (W-W_copy).norm();
    std::cout << "iteration = " << i+1 << ", rel. residual = " << res << ", maximum error = " << Err.col(index).norm() << ", solution incr. = " << inc << std::endl;
    if(inc < tol_) break;
    W_copy = W;
*/
  }

  double Wnorm = W.froNorm();
  for(int i=0; i<communicator_->numCPUs(); ++i) { if(i==communicator_->myID()) {
  std::cerr << "norm of W = " << Wnorm << std::endl;
  } communicator_->sync(); }
  
  for(int j=1; j<=rowCount_; j++) {
    int p = _FORTRAN(indxg2p)(&j, &blockSize_, &dummy, &zero, &rowCpus);
    if(myrow == p) {
      int lj = _FORTRAN(indxg2l)(&j, &blockSize_, NULL, NULL, &rowCpus) - 1;
      W.getMatrixRow(j, row);
      for(int k=0; k<basisDimension_; ++k) basisBuffer_(lj,k) = row[k];
    }
  }

  delete [] row;

  Cblacs_gridexit(context);
  Cfree_blacs_system_handle(blacsHandle);
#endif
}

void
DistrNonnegativeMatrixFactorization::solveNNLS_MRHS(SCDoubleMatrix &A, SCDoubleMatrix &B, SCDoubleMatrix &X, int flag)
{
  // if flag = 0 then solve: min ||AX^T-B|| s.t. X >= 0
  // if flag = 1 then solve: min ||AX^T-B^T|| s.t. X >= 0

  Plh solver(A.getNumberOfRows(), A.getNumberOfCols());
  solver.init();
  solver.setMatrix(A);

  double *rhs = new double[A.getNumberOfRows()];
  double *sol = new double[A.getNumberOfCols()];
  int nrhs = X.getNumberOfRows();

  for(int i=1; i<=nrhs; ++i) {
    if(flag==0) B.getMatrixColumn(i,rhs);
    else        B.getMatrixRow(i,rhs);
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
