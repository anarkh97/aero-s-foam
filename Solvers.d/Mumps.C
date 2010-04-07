#include <Driver.d/Communicator.h>
#include <map>

#ifdef USE_MUMPS
#include <dmumps_c.h> // double precision mumps header
#include <zmumps_c.h> // complex double precision 

inline void Tmumps_c(DMUMPS_STRUC_C &id) { dmumps_c(&id); }
inline void Tmumps_c(ZMUMPS_STRUC_C &id) { zmumps_c(&id); }
#endif

#include <Driver.d/Domain.h>
extern Domain * domain;

#define	USE_COMM_WORLD	-987654 
// MUMPS coded in Fortran -> watch out for index to match documentation
#define ICNTL(I) 	icntl[(I)-1]
#define CNTL(I)        	cntl[(I)-1]
#define INFOG(I)        infog[(I)-1]
#define INFO(I)         info[(I)-1]

template<class Scalar>
GenMumpsSolver<Scalar>::GenMumpsSolver(Connectivity *nToN, EqNumberer *_dsa, int *map, FSCommunicator *_mpicomm)
 : SparseData(_dsa,nToN,map,0,1)
{
  neq = numUncon;
  nNonZero = xunonz[numUncon]-1; 
  //cerr << "int Mumps constructor #1, mumps_sym = " << domain->solInfo().mumps_sym << ", neq = " << neq << ", nNonZero = " << nNonZero << endl;
  unonz = new Scalar[nNonZero];
  for(int i = 0; i < nNonZero; ++i) unonz[i] = 0.0;
  mpicomm = _mpicomm;
  init();
}

template<class Scalar>
GenMumpsSolver<Scalar>::GenMumpsSolver(Connectivity *nToN, DofSetArray *_dsa, ConstrainedDSA *c_dsa, FSCommunicator *_mpicomm)
 : SparseData(_dsa,c_dsa,nToN,0,1)
{
  neq = numUncon;
  myMem = 0; 
  nNonZero = xunonz[numUncon]-1;
  //cerr << "int Mumps constructor #2, mumps_sym = " << domain->solInfo().mumps_sym << ", neq = " << neq << ", nNonZero = " << nNonZero << endl;
  unonz	= new Scalar[nNonZero];
  for(int i = 0; i < nNonZero; ++i) unonz[i] = 0.0;
  mpicomm = _mpicomm;
  init();
}

template<class Scalar>
void
GenMumpsSolver<Scalar>::init()
{
#ifdef USE_MUMPS
  mumpsId.id.par = 1; // 1: working host model
  mumpsId.id.sym = domain->solInfo().pivot ? 2 : 1; // 2: general symmetric, 1: symmetric positive definite, 0: unsymmetric 
#ifdef USE_MPI
  //mumpsId.id.comm_fortran = USE_COMM_WORLD; // default value for fortran communicator
  if(mpicomm) mumpsId.id.comm_fortran = MPI_Comm_c2f(mpicomm->getComm());
  else mumpsId.id.comm_fortran = MPI_Comm_c2f(MPI_COMM_SELF);
#else
  mumpsId.id.comm_fortran = USE_COMM_WORLD; // default value for fortran communicator
#endif
  mumpsId.id.job = -1; // initialize instance of the mumps package
  Tmumps_c(mumpsId.id);
  host = (mpicomm) ? (mpicomm->cpuNum() == 0) : true;

  // Set control parameters CNTL and ICNTL
  map<int,double>::iterator CntlIter = domain->solInfo().mumps_cntl.begin();
  while(CntlIter != domain->solInfo().mumps_cntl.end()) {
    int CntlNum         = CntlIter->first;
    double CntlPar      = CntlIter->second;
    mumpsId.id.CNTL(CntlNum) = CntlPar;
    CntlIter ++;
  }
  map<int,int>::iterator IcntlIter = domain->solInfo().mumps_icntl.begin();
  while(IcntlIter != domain->solInfo().mumps_icntl.end()) {
    int IcntlNum        = IcntlIter->first;
    int IcntlPar        = IcntlIter->second;
    mumpsId.id.ICNTL(IcntlNum) = IcntlPar;
    IcntlIter ++;
  }
  // NOTE: centralized assembled matrix input (ICNTL(5) = 0 with ICNTL(18) = 0) is the only option fully supported
  //       distributed assembled matrix input (ICNTL(5) = 0 with ICNTL(18) = 3) is currently supported only for the FETI-DP coarse_solver
  if(mumpsId.id.ICNTL(5) != 0) {
    cerr << "user defined ICNTL(5) not supported, setting to 0\n";
    mumpsId.id.ICNTL(5) = 0; // 0: assembled matrix input
  }
  if(mumpsId.id.ICNTL(18) != 0 && mumpsId.id.ICNTL(18) != 3) {
    cerr << "user defined ICNTL(18) not supported, setting to 0\n";
    mumpsId.id.ICNTL(18) = 0; // 0: centralized assembled matrix input, 3: distributed assembled matrix input
  }
  // NOTE: centralized dense right-hand-side and solution (ICNTL(20) = 0 with ICNTL(21) = 0) is the only option supported
  if(mumpsId.id.ICNTL(20) != 0) {
    cerr << "user defined ICNTL(20) not supported, setting to 0\n";
    mumpsId.id.ICNTL(20) = 0; // 0: dense RHS 
  }
  if(mumpsId.id.ICNTL(21) != 0) {
    cerr << "user defined ICNTL(21) not supported, setting to 0\n";
    mumpsId.id.ICNTL(21) = 0; // 0: centralized solution
  }
  if(domain->solInfo().pivot) { // matrix is not assumed to be positive definite, may be singularities 
    mumpsId.id.ICNTL(24) = 1; // 1: enable null pivot row detection
    mumpsId.id.ICNTL(13) = 1; // 1: ScaLAPACK will not be used for the root frontal matrix (recommended for null pivot row detection)
    mumpsId.id.CNTL(3) = -domain->solInfo().trbm; // tolerance used to detect zero pivots during factorization (not used for SPD matrix)
  }
#endif
}

template<class Scalar>
void
GenMumpsSolver<Scalar>::zeroAll()
{
  if(!unonz) unonz = new Scalar[nNonZero];
  for(int i = 0; i < nNonZero; ++i) unonz[i] = 0.0;
}

template<class Scalar>
void
GenMumpsSolver<Scalar>::add(FullSquareMatrix &kel, int *dofs)
{
  int i, j, m, mstart, mstop;
  int kndof = kel.dim();                       // Dimension of element stiff.

  for(i = 0; i < kndof; ++i) {                 // Loop over rows.
    if(unconstrNum[dofs[i]] == -1) continue;   // Skip constrained dofs
    for(j = 0; j < kndof; ++j) {               // Loop over columns.
      if(unconstrNum[dofs[j]] == -1) continue; // Skip constrained dofs
      if(unconstrNum[dofs[j]] < unconstrNum[dofs[i]]) continue;
      mstart = xunonz[unconstrNum[dofs[j]]];
      mstop  = xunonz[unconstrNum[dofs[j]]+1];
      for(m = mstart; m < mstop; ++m) {
        if(rowu[m-1] == (unconstrNum[dofs[i]] + 1)) {
          unonz[m-1] += kel[i][j];
          break;
        }
      }
    }
  }
}

template<class Scalar>
void
GenMumpsSolver<Scalar>::add( int dofi, int dofj, Scalar d)
{
 // WARNING: this adds only [i,j] with  j >= i (i.e. upper part)
  if((dofi < 0) || (dofj < 0)) return;
  int m, mstart, mstop, rowi, colj;
  if(unconstrNum) {
    if((rowi = unconstrNum[dofi]) == -1 || (colj = unconstrNum[dofj]) == -1) return;
  }
  else { rowi = dofi; colj = dofj; }

  if(colj<rowi) { // swap row & col to be in the upper part
    int tmp = colj;
    colj = rowi;
    rowi = tmp;
  }
  // upper part
  mstart = xunonz[colj];
  mstop  = xunonz[colj+1];
  for(m = mstart; m < mstop; ++m) {
    if(rowu[m-1] == (rowi+1)) {
      unonz[m-1] += d;
      break;
    }
  }
}

template<class Scalar>
void
GenMumpsSolver<Scalar>::add(GenAssembledFullM<Scalar> &kel, int *dofs)
{
  // this function is used to assemble Kcc and requires dofs to be in constrained numbering
  int i, j, m, mstart, mstop, ri, rj;
  for(i = 0; i < kel.numRow(); ++i) {       // Loop over rows.
    if((ri = dofs[i]) == -1) continue;      // Skip constrained dofs
    for(j = 0; j < kel.numCol(); ++j) {     // Loop over columns.
      if((rj = dofs[j]) == -1) continue;    // Skip constrained dofs
      if(rj < ri) continue;
      mstart = xunonz[rj];
      mstop  = xunonz[rj+1];
      for(m = mstart; m < mstop; ++m) {
        if(rowu[m-1] == (ri + 1)) {
          unonz[m-1] += kel[i][j];
          break;
        }
      }
    }
  }
}

template<class Scalar>
void
GenMumpsSolver<Scalar>::add(GenFullM<Scalar> &kel, int fi, int fj)
{
  int i, j, m, mstart, mstop;
  for(i = 0; i < kel.numRow(); ++i ) {      // Loop over rows.
    if(unconstrNum[fi+i] == -1) continue;   // Skip constrained dofs
    for(j = 0; j < kel.numCol(); ++j) {     // Loop over columns.
      if(unconstrNum[fj+j] == -1) continue; // Skip constrained dofs
      if(unconstrNum[fj+j] < unconstrNum[fi+i]) continue;
      mstart = xunonz[unconstrNum[fj+j]];
      mstop  = xunonz[unconstrNum[fj+j]+1];
      for(m = mstart; m < mstop; ++m) {
        if(rowu[m-1] == (unconstrNum[fi+i] + 1) ) {
          unonz[m-1] += kel[i][j];
          break;
        }
      }
    }
  }
}

template<class Scalar>
Scalar
GenMumpsSolver<Scalar>::diag(int dof) const
{
  int m, mstart, mstop;

  mstart = xunonz[dof]-1;
  mstop  = xunonz[dof+1]-1;

  for(m = mstart; m < mstop; ++m) {
    if(rowu[m]-1 == dof) {
      if(unonz[m] == 0.0)
        return (1.0);
      else
        return unonz[m];
    }
  }
  throw "GenMumpsSolver<Scalar>::diag - 1 - this should never be reached";
}

template<class Scalar>
Scalar &
GenMumpsSolver<Scalar>::diag(int dof)
{
  int m, mstart, mstop;

  mstart = xunonz[dof]-1;
  mstop  = xunonz[dof+1]-1;

  for(m = mstart; m < mstop; ++m) {
    if(rowu[m]-1 == dof) {
        return unonz[m];
    }
  }
  throw "GenMumpsSolver<Scalar>::diag - 2 - this should never be reached";
}

template<class Scalar>
void
GenMumpsSolver<Scalar>::addDiscreteMass(int dof, Scalar dmass)
{
  if(dof < 0) return;
  int cdof;
  if(unconstrNum) cdof = unconstrNum[dof]; 
  else cdof = dof;
  if(cdof < 0) return;

  int mstart = xunonz[cdof];
  int mstop  = xunonz[cdof+1];
  for(int m = mstart; m<mstop; ++m) {
    if(rowu[m-1] == (cdof + 1)) {
      unonz[m-1] += dmass;
      break;
    }
  }
}

template<class Scalar>
void
GenMumpsSolver<Scalar>::unify(FSCommunicator *_communicator)
{
#if defined(USE_MPI) && defined(USE_MUMPS)
  if(mpicomm && mumpsId.id.ICNTL(18) == 0)
    mpicomm->reduce(xunonz[numUncon]-1, unonz, 0); // assemble on host (mpi process with rank 0)
#endif
}

template<class Scalar>
void
GenMumpsSolver<Scalar>::factor()
{
  if(numUncon == 0) return;
#ifdef USE_MUMPS
  if(host) {
    mumpsId.id.n = neq;
    if(mumpsId.id.ICNTL(18) == 0) { // centralized matrix input
      mumpsId.id.nz  = nNonZero;
      mumpsId.id.irn = rowu; //irn;
      mumpsId.id.jcn = colu; //jcn;
      copyToMumpsLHS(mumpsId.id.a, unonz, nNonZero);
    }
  }
  if(mumpsId.id.ICNTL(18) == 3) { // distributed matrix input
    mumpsId.id.nz_loc  = nNonZero;
    mumpsId.id.irn_loc = rowu; //irn;
    mumpsId.id.jcn_loc = colu; //jcn;
    copyToMumpsLHS(mumpsId.id.a_loc, unonz, nNonZero);
  }
  mumpsId.id.job = 4; // 4: analysis and factorization, 1: analysis only, 2: factorization only
  Tmumps_c(mumpsId.id);

  if(host && mumpsId.id.INFOG(1) < 0) {
    cerr << "error in MUMPS factorization: INFOG(1) = " << mumpsId.id.INFOG(1) << ", INFOG(2) = " << mumpsId.id.INFOG(2) << endl;
    exit(-1);
  }

  nrbm = mumpsId.id.INFOG(28); // number of zero pivots detected
  if(this->print_nullity && host && nrbm > 0) 
    cerr << " ... Matrix is singular: size = " << neq << ", rank = " << neq-nrbm << ", nullity = " << nrbm << " ...\n";

#endif
}

template<class Scalar>
void
GenMumpsSolver<Scalar>::solve(Scalar *rhs, Scalar *solution)
{
  this->solveTime = -getTime();
#ifdef USE_MUMPS
  if(host) {
    copyToMumpsRHS(mumpsId.id.rhs, rhs, numUncon); // mumpsId.id.rhs = copy of rhs;
    mumpsId.id.nrhs = 1;
  }
  mumpsId.id.job = 3; // 3: solve
  Tmumps_c(mumpsId.id);
  if(host) copyFromMumpsRHS(solution, mumpsId.id.rhs, numUncon); // solution = mumpsId.id.rhs;
#ifdef USE_MPI
  if(mpicomm) mpicomm->broadcast(numUncon, solution, 0); // send from host to others
#endif
#endif
  this->solveTime += getTime();
}

template<class Scalar>
void
GenMumpsSolver<Scalar>::reSolve(int nRHS, Scalar *rhs)
{
  this->solveTime = -getTime();
#ifdef USE_MUMPS
  if(host) {
    copyToMumpsRHS(mumpsId.id.rhs, rhs, numUncon*nRHS); // mumpsId.id.rhs = copy of rhs;
    mumpsId.id.nrhs = nRHS;
    mumpsId.id.lrhs = numUncon; // leading dimension 
  }
  mumpsId.id.job = 3; // 3: solve
  Tmumps_c(mumpsId.id);
  if(host) copyFromMumpsRHS(rhs, mumpsId.id.rhs, numUncon*nRHS); // rhs = mumpsId.id.rhs;
#ifdef USE_MPI
  if(mpicomm) mpicomm->broadcast(numUncon*nRHS, rhs, 0); // send from host to all others
#endif
#endif
  this->solveTime += getTime();
}

template<class Scalar>
void
GenMumpsSolver<Scalar>::reSolve(int nRHS, Scalar **rhs)
{
  this->solveTime = -getTime();
#ifdef USE_MUMPS
  if(host) {
    copyToMumpsRHS(mumpsId.id.rhs, rhs, numUncon, nRHS); // mumpsId.id.rhs = copy of rhs;
    mumpsId.id.nrhs = nRHS;
    mumpsId.id.lrhs = numUncon; // leading dimension 
  }
  mumpsId.id.job = 3; // 3: solve
  Tmumps_c(mumpsId.id);
  if(host) copyFromMumpsRHS(rhs, mumpsId.id.rhs, numUncon, nRHS); // rhs = mumpsId.id.rhs;
#ifdef USE_MPI
  if(mpicomm) for(int i=0; i<nRHS; ++i) mpicomm->broadcast(numUncon, rhs[i], 0); // send from host to others
#endif
#endif
  this->solveTime += getTime();
}

template<class Scalar>
void
GenMumpsSolver<Scalar>::reSolve(int nRHS, GenVector<Scalar> *rhs)
{
  this->solveTime = -getTime();
#ifdef USE_MUMPS
  if(host) {
    copyToMumpsRHS(mumpsId.id.rhs, rhs, numUncon, nRHS); // mumpsId.id.rhs = copy of rhs;
    mumpsId.id.nrhs = nRHS;
    mumpsId.id.lrhs = numUncon; // leading dimension 
  }
  mumpsId.id.job = 3; // 3: solve
  Tmumps_c(mumpsId.id);
  if(host) copyFromMumpsRHS(rhs, mumpsId.id.rhs, numUncon, nRHS); // rhs = mumpsId.id.rhs;
#ifdef USE_MPI
  if(mpicomm) for(int i=0; i<nRHS; ++i) mpicomm->broadcast(numUncon, rhs[i].data(), 0); // send from host to others
#endif
#endif
  this->solveTime += getTime();
}

template<class Scalar>
void
GenMumpsSolver<Scalar>::getNullSpace(Scalar *rbms)
{
#ifdef USE_MUMPS
  mumpsId.id.ICNTL(25) = -1;
  reSolve(nrbm, rbms);
  mumpsId.id.ICNTL(25) = 0;
#endif
}

template<class Scalar>
void
GenMumpsSolver<Scalar>::print()
{
#ifdef USE_MUMPS
 if(host && mumpsId.id.ICNTL(18) == 0) { // centralized matrix on host
   cerr << endl;
   for(int i=0; i < nNonZero; ++i)
     cerr << "K(" << mumpsId.id.irn[i] << "," << mumpsId.id.jcn[i] << ") = " << unonz[i] << endl;
 }
 else if(host && mumpsId.id.ICNTL(18) == 3) cerr << "GenMumpsSolver<Scalar>::print() is not implemented for ICNTL(18) = 3\n";
#endif
}

template<class Scalar>
long int
GenMumpsSolver<Scalar>::size()
{
#ifdef USE_MUMPS
#ifdef DISTRIBUTED
  if(mumpsId.id.ICNTL(18) == 3) 
    return mpicomm->globalSum(nNonZero);
  else 
#endif
#endif
  return nNonZero;
}

template<class Scalar>
GenMumpsSolver<Scalar>::~GenMumpsSolver()
{
  if(unonz) { delete [] unonz; unonz = 0; }
#ifdef USE_MUMPS
  mumpsId.id.job = -2; // -2: destroys instance of mumps package
  Tmumps_c(mumpsId.id);
#endif
}

template<>
void GenMumpsSolver<complex<double> >::addImaginary(FullSquareMatrix &kel, int *dofs);

template<>
void GenMumpsSolver<double>::addImaginary(FullSquareMatrix &kel, int *dofs);

#ifdef USE_MUMPS
template<>
void
GenMumpsSolver<DComplex>::copyToMumpsRHS(mumps_double_complex *&m, DComplex *d, int len);

template<>
void
GenMumpsSolver<DComplex>::copyFromMumpsRHS(DComplex *d, mumps_double_complex *m, int len);

template<>
void
GenMumpsSolver<double>::copyToMumpsRHS(double *&m, double *d, int len);

template<>
void
GenMumpsSolver<double>::copyFromMumpsRHS(double *d, double *m, int len);

template<>
void
GenMumpsSolver<DComplex>::copyToMumpsLHS(mumps_double_complex *&m, DComplex *&d, int len);

template<>
void
GenMumpsSolver<double>::copyToMumpsLHS(double *&m, double *&d, int len);
#endif
