#include <Driver.d/Communicator.h>
#include <map>

#ifdef USE_MUMPS

inline void Tmumps_c(DMUMPS_STRUC_C &id) { dmumps_c(&id); }
inline void Tmumps_c(ZMUMPS_STRUC_C &id) { zmumps_c(&id); }
#endif

#include <Driver.d/Domain.h>
extern Domain * domain;
extern long totMemMumps;

#define	USE_COMM_WORLD	-987654 
// MUMPS coded in Fortran -> watch out for index to match documentation
#define ICNTL(I) 	icntl[(I)-1]
#define CNTL(I)        	cntl[(I)-1]
#define INFOG(I)        infog[(I)-1]
#define INFO(I)         info[(I)-1]
#define RINFOG(I)       rinfog[(I)-1]
#define RINFO(I)        rinfo[(I)-1]

template<class Scalar>
GenMumpsSolver<Scalar>::GenMumpsSolver(Connectivity *nToN, EqNumberer *_dsa, int *map, FSCommunicator *_mpicomm)
 : SparseData(_dsa,nToN,map,0,1)
{
#ifndef USE_MUMPS
  std::cerr << " *** ERROR: Solver requires AERO-S configured with the MUMPS library. Exiting...\n";
  exit(-1);
#endif
  neq = numUncon;
  nNonZero = xunonz[numUncon]-1; 
  unonz = new Scalar[nNonZero];
  for(int i = 0; i < nNonZero; ++i) unonz[i] = 0.0;
  mpicomm = _mpicomm;
  init();
}

template<class Scalar>
GenMumpsSolver<Scalar>::GenMumpsSolver(Connectivity *nToN, DofSetArray *_dsa, ConstrainedDSA *c_dsa, FSCommunicator *_mpicomm)
 : SparseData(_dsa,c_dsa,nToN,0,1,domain->solInfo().unsym())
{
#ifndef USE_MUMPS
  std::cerr << " *** ERROR: Solver requires AERO-S configured with the MUMPS library. Exiting...\n";
  exit(-1);
#endif
  neq = numUncon;
  myMem = 0; 
  nNonZero = xunonz[numUncon]-1;
  unonz	= new Scalar[nNonZero];
  for(int i = 0; i < nNonZero; ++i) unonz[i] = 0.0;
  mpicomm = _mpicomm;
  init();
}

template<class Scalar>
GenMumpsSolver<Scalar>::GenMumpsSolver(Connectivity *nToN, DofSetArray *_dsa, ConstrainedDSA *c_dsa, int nsub,
                                       GenSubDomain<Scalar> **sd, FSCommunicator *_mpicomm)
 : SparseData(_dsa,c_dsa,nToN,0,1,domain->solInfo().unsym()), MultiDomainSolver<Scalar>(numUncon, nsub, sd, _mpicomm)
{
#ifndef USE_MUMPS
  std::cerr << " *** ERROR: Solver requires AERO-S configured with the MUMPS library. Exiting...\n";
  exit(-1);
#endif
  neq = numUncon;
  myMem = 0;
  nNonZero = xunonz[numUncon]-1;
  unonz = new Scalar[nNonZero];
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
  if(domain->solInfo().unsym()) mumpsId.id.sym = 0;
#ifdef USE_MPI
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
  }
#endif
}

template<class Scalar>
void
GenMumpsSolver<Scalar>::zeroAll()
{
  if(!unonz) unonz = new Scalar[nNonZero];
  for(int i = 0; i < nNonZero; ++i) unonz[i] = 0.0;
#ifdef USE_MUMPS
  mumpsId.id.job = -2; // destroy instance of the mumps package
  Tmumps_c(mumpsId.id);
#endif
  init();
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
      if(!domain->solInfo().unsym() && unconstrNum[dofs[j]] < unconstrNum[dofs[i]]) continue;
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
      mumpsId.id.irn = rowu;
      mumpsId.id.jcn = colu;
      copyToMumpsLHS(mumpsId.id.a, unonz, nNonZero);
    }
  }
  if(mumpsId.id.ICNTL(18) == 3) { // distributed matrix input
    mumpsId.id.nz_loc  = nNonZero;
    mumpsId.id.irn_loc = rowu;
    mumpsId.id.jcn_loc = colu;
    copyToMumpsLHS(mumpsId.id.a_loc, unonz, nNonZero);
  }
  mumpsId.id.job = 4; // 4: analysis and factorization, 1: analysis only, 2: factorization only
  while(true) {
    Tmumps_c(mumpsId.id);

    if(mumpsId.id.INFOG(1) == -8 || mumpsId.id.INFOG(1)  == -9) { 
       if(host) cerr << " ... increasing MUMPS workspace     ...\n"; 
       mumpsId.id.ICNTL(14) *= 2; 
       mumpsId.id.job = 2; // recall factorization
    }
    else if(mumpsId.id.INFOG(1) < 0) { 
      if(host) cerr << " *** ERROR: MUMPS factorization returned error code. Exiting...\n";
      exit(-1);
    }
    else break;
  }

  nrbm = mumpsId.id.INFOG(28); // number of zero pivots detected
  if(this->print_nullity && host && nrbm > 0) 
    cerr << " ... Matrix is singular: size = " << neq << ", rank = " << neq-nrbm << ", nullity = " << nrbm << " ...\n";

  totMemMumps += mumpsId.id.INFOG(19)*1024; // INFOG(19) is the size in millions of bytes of all mumps internal data 
                                            // allocated during factorization: sum over all processors
#endif
}

template<class Scalar>
int*
GenMumpsSolver<Scalar>::getPivnull_list()
{
#ifdef USE_MUMPS
  return mumpsId.id.pivnul_list;
#else
  return NULL;
#endif
} 

template<class Scalar>
void
GenMumpsSolver<Scalar>::solve(Scalar *rhs, Scalar *solution)
{
  if(numUncon == 0) return;
  this->solveTime -= getTime();
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
  if(mumpsId.id.ICNTL(11) > 0) printStatistics();
#endif
  this->solveTime += getTime();
}

template<class Scalar>
void
GenMumpsSolver<Scalar>::reSolve(int nRHS, Scalar *rhs)
{
  if(numUncon == 0) return;
  this->solveTime -= getTime();
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
  if(mumpsId.id.ICNTL(11) > 0) printStatistics();
#endif
  this->solveTime += getTime();
}

template<class Scalar>
void
GenMumpsSolver<Scalar>::reSolve(int nRHS, Scalar **rhs)
{
  if(numUncon == 0) return;
  this->solveTime -= getTime();
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
  if(mumpsId.id.ICNTL(11) > 0) printStatistics();
#endif
  this->solveTime += getTime();
}

template<class Scalar>
void
GenMumpsSolver<Scalar>::reSolve(int nRHS, GenVector<Scalar> *rhs)
{
  if(numUncon == 0) return;
  this->solveTime -= getTime();
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
  if(mumpsId.id.ICNTL(11) > 0) printStatistics();
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
 if((mumpsId.id.ICNTL(18) == 0 && host) || mumpsId.id.ICNTL(18) == 3) { // centralized matrix on host
   for(int i = 0; i < nNonZero; ++i)
     cerr << "A(" << rowu[i] << "," << colu[i] << ") = " << unonz[i] << endl;
 }
#endif
}

template<class Scalar>
long int
GenMumpsSolver<Scalar>::size()
{
#ifdef USE_MUMPS
  if(mumpsId.id.INFOG(3) < 0) return static_cast<long int>(-1e6)*mumpsId.id.INFOG(3);
  else return mumpsId.id.INFOG(3);
#else
  return 0;
#endif
}

template<class Scalar>
void
GenMumpsSolver<Scalar>::printStatistics()
{
#ifdef USE_MUMPS
  if(host) {
    cerr << "RINFOG(4,...,11) = " << "("
         << mumpsId.id.RINFOG(4) << ", " << mumpsId.id.RINFOG(5) << ", "  << mumpsId.id.RINFOG(6) << ", "
         << mumpsId.id.RINFOG(7) << ", " << mumpsId.id.RINFOG(8) << ", "  << mumpsId.id.RINFOG(9) << ", "
         << mumpsId.id.RINFOG(10) << ", " << mumpsId.id.RINFOG(11) << ")" << endl;
  }
#endif
}

template<class Scalar>
GenMumpsSolver<Scalar>::~GenMumpsSolver()
{
#ifdef USE_MUMPS
  if(mpicomm) mpicomm->sync();

  if(host) {
    if(mumpsId.id.ICNTL(18) == 0) { // centralized matrix input
      if((void*)mumpsId.id.a != (void*)unonz) {
        if (mumpsId.id.a) delete[] mumpsId.id.a;
        mumpsId.id.a = 0;
      }
    }
  }

  if(mumpsId.id.ICNTL(18) == 3) {
    if((void*)mumpsId.id.a_loc != (void*)unonz) {
      if (mumpsId.id.a_loc) delete[] mumpsId.id.a_loc;
      mumpsId.id.a_loc = 0;
    }
  }

  if(unonz) {
    delete [] unonz;
    unonz = 0;
   }

  mumpsId.id.job = -2; // -2: destroys instance of mumps package
  Tmumps_c(mumpsId.id);

#endif
  if(mpicomm) mpicomm->sync();
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
