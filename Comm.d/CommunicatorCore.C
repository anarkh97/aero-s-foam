#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <Utils.d/linkfc.h>
#include <Utils.d/dbg_alloca.h>
#include <Comm.d/Communicator.h>

#ifdef TFLOP
#include "mpi.h"
#define SEGMENT_GLBSUM
#endif

FILE *debugFile = 0;

#ifdef USE_MPI
#ifndef MPI_INTEGER
#define MPI_INTEGER MPI_INT
#endif

#ifndef MPI_CHARACTER
#define MPI_CHARACTER MPI_CHAR
#endif

// PJSA: for linux
#if (defined(USE_MPI) && !defined(MPI_DOUBLE_COMPLEX) && defined(LAM_MPI))
//#include<mpisys.h>
#define MPI_DOUBLE_COMPLEX ((MPI_Datatype) &lam_mpi_cxx_dblcplex)
#endif
// CRW - this define is a hack to let it compile on sun, but feti-h will not work
#if (defined(USE_MPI) && !defined(MPI_DOUBLE_COMPLEX) && defined(NO_COMPLEX_MPI))    //CRW
#define MPI_DOUBLE_COMPLEX MPI_LONG_DOUBLE    //CRW
#endif    //CRW

#ifdef CXX  // for DEC cxx compiler
MPI_Datatype CommTrace<double>::MPIType = MPI_DOUBLE;
MPI_Datatype CommTrace<float>::MPIType = MPI_FLOAT;
MPI_Datatype CommTrace<int>::MPIType = MPI_INTEGER;
MPI_Datatype CommTrace<char>::MPIType = MPI_CHARACTER;
MPI_Datatype CommTrace<DComplex>::MPIType = MPI_DOUBLE_COMPLEX;
//MPI_Datatype CommTrace<bool>::MPIType = MPI_BOOL;
#else 
template<>
MPI_Datatype CommTrace<double>::MPIType = MPI_DOUBLE;
template<>
MPI_Datatype CommTrace<float>::MPIType = MPI_FLOAT;
template<>
MPI_Datatype CommTrace<int>::MPIType = MPI_INTEGER;
template<>
MPI_Datatype CommTrace<char>::MPIType = MPI_CHARACTER;
template<>
MPI_Datatype CommTrace<DComplex>::MPIType = MPI_DOUBLE_COMPLEX;
//template<>
//MPI_Datatype CommTrace<bool>::MPIType = MPI_BOOL;
#endif
#endif

static MPI_Request nullReq;
static MPI_Status nullStat;

//------------------------------------------------------------------------------ 

Communicator::Communicator(int _ncpu) 
: pendReq(nullReq), reqStatus(nullStat), nPendReq(0)
{
  comm = MPI_COMM_WORLD;
  nPendReq = 0;
  glNumCPU = _ncpu;
}

void
Communicator::sync()
{
#ifdef USE_MPI
  MPI_Barrier(comm);
#endif
}

int
Communicator::myID()
{
  int id = 0;
#ifdef USE_MPI
  MPI_Comm_rank(comm, &id);
#endif
  return id;
}

SysCom::SysCom(int &argc, char **&argv) 
: Communicator(MPI_COMM_WORLD, stderr)
{
  salinasCommunicator = 0;
#if defined(USE_MPI) && !defined(CREATE_DSO)
  MPI_Init(&argc,&argv);
#endif
}

SysCom::SysCom(MPI_Comm *mpi_comm) : Communicator(*mpi_comm, stderr)
{
  salinasCommunicator = mpi_comm;
#ifdef TFLOP
  int structID = 2;
  _FORTRAN(startcom)(mpi_comm, structID);
#endif
}

SysCom::SysCom() 
: Communicator(MPI_COMM_WORLD,stderr)
{
  salinasCommunicator = 0;
}

SysCom::~SysCom()
{
#if defined(USE_MPI) && !defined(CREATE_DSO)
  if(!salinasCommunicator) MPI_Finalize();
#endif
}

int
Communicator::numCPUs()
{
  int id = 1;
#ifdef USE_MPI
  MPI_Comm_size(comm, &id);
#endif
  return id;
}

int CPairCompare(const void *a, const void *b)
{
 CPair *pa = (CPair *)a;
 CPair *pb = (CPair *)b;

 int mina, minb, maxa, maxb;
 if(pa->glFrom < pa->glTo) {
   mina = pa->glFrom;
   maxa =  pa->glTo;
 } else {
  mina = pa->glTo;
  maxa = pa->glFrom;
 }

 if(pb->glFrom < pb->glTo) {
   minb = pb->glFrom;
   maxb =  pb->glTo;
 } else {
  minb = pb->glTo;
  maxb = pb->glFrom;
 }

 if(mina < minb || (mina == minb && maxa < maxb)) return -1;
 if(mina == minb && maxa == maxb) return 0;
 return 1;
}

static int v = 0;

int uniqueTag()
{
 v = (v+1)%100;
 return v+100;
}

void Communicator::waitForAllReq()
{
#ifdef USE_MPI
  // Just making sure that reqStatus has an appropriate length
  MPI_Status *safe = reqStatus+nPendReq;
 
  if (safe == 0)
    exit(1);
 
  int nSuccess = MPI_Waitall(nPendReq, pendReq+0, reqStatus+0);
 
  if (nSuccess) {
    fprintf(stderr, " *** ERROR: unexpected success number %d\n", nSuccess);
    exit(1);
  }
 
  nPendReq = 0;
#endif
}  

void Communicator::split(int color, int maxcolor, Communicator** c)
{
  int i;
  for (i=0; i<maxcolor; ++i)
    c[i] = 0;

#ifdef USE_MPI
  int rank;
  MPI_Comm_rank(comm, &rank);
  MPI_Comm comm1;
  MPI_Comm_split(comm, color+1, rank, &comm1); //wrong should use color
  c[color] = new Communicator(comm1,stderr);

  int* leaders = new int[maxcolor];
  int* newleaders = new int[maxcolor];
  for (i=0; i<maxcolor; ++i)
    leaders[i] = -1;

  int localRank;
  MPI_Comm_rank(comm1, &localRank);
  if (localRank == 0)
    leaders[color] = rank;

  MPI_Allreduce(leaders, newleaders, maxcolor, MPI_INTEGER, MPI_MAX, comm);

   for (i=0; i<maxcolor; ++i) {
     if (i != color && newleaders[i] >= 0) {
       int tag;
       if (color < i)
        tag = maxcolor*(color+1)+i+1;
       else
        tag = maxcolor*(i+1)+color+1;
       MPI_Comm comm2;
       MPI_Intercomm_create(comm1, 0, comm, newleaders[i], tag, &comm2);
       c[i] = new Communicator(comm2,stderr);
     }
   }

  if (leaders) 
    delete [] leaders;
  if (newleaders) 
    delete [] newleaders;
#else
  c[color] = this;
#endif
}

Communicator::Communicator(MPI_Comm c1, FILE *fp)
  : pendReq(nullReq), reqStatus(nullStat)
{
  comm = c1;
  nPendReq = 0;
}

Communicator::Communicator(const Communicator &c1)
: pendReq(nullReq) , reqStatus(nullStat), nPendReq(0)
{
  comm = c1.comm;
#ifdef USE_MPI
  MPI_Comm_size(comm, &glNumCPU);
#else
  glNumCPU = 1;
#endif
}

int
Communicator::remoteSize()
{
  int numRemote = 0;
#ifdef USE_MPI 
  MPI_Comm_remote_size(comm, &numRemote);
#endif
  return numRemote;
}

bool Communicator::globalMax(bool b)
{
  int buff;
#ifdef USE_MPI
  int data = (b) ? 1 : 0;
  MPI_Allreduce(&data, &buff, 1, CommTrace<int>::MPIType, MPI_MAX, comm);
#endif
  return buff != 0;
}

#ifdef NO_COMPLEX_MPI
template <>
complex<double> Communicator::globalSum(complex<double> data)
{
  double tmp_data[2] = { data.real(), data.imag() };
  globalSum(2, tmp_data);
  return complex<double>(tmp_data[0],tmp_data[1]);
}

template <>
complex<double> Communicator::globalMax(complex<double> data)
{
  cerr << "ERROR: Communicator::globalMax called with complex data\n";
  return complex<double>(0.0,0.0);
}

template <>
complex<double> Communicator::globalMin(complex<double> data)
{
  cerr << "ERROR: Communicator::globalMin called with complex data\n";
  return complex<double>(0.0,0.0);
}

template <>
void Communicator::globalSum(int num, complex<double> *data)
{
  double *tmp_data = new double[2*num];
  for(int i=0; i<num; ++i) { tmp_data[2*i] = data[i].real(); tmp_data[2*i+1] = data[i].imag(); }
  globalSum(2*num,tmp_data);
  for(int i=0; i<num; ++i) data[i] = complex<double>(tmp_data[2*i],tmp_data[2*i+1]);
  delete [] tmp_data;
}

template <>
void Communicator::sendTo(int cpu, int tag, complex<double> *buffer, int len)
{
  cerr << "ERROR: Communicator::sendTo called with complex data\n";
}
  
template <>
RecInfo Communicator::recFrom(int tag, complex<double> *buffer, int len)
{
  cerr << "ERROR: Communicator::recFrom called with complex data\n";
  return RecInfo();
}
  
template <>
void Communicator::allGather(complex<double> *send_data, int send_count, complex<double> *recv_data, int recv_count)
{
  cerr << "ERROR: Communicator::allGather called with complex data\n";
}

template <>
void Communicator::allGatherv(complex<double> *send_data, int send_count, complex<double> *recv_data, int recv_counts[], int displacements[])
{ 
  cerr << "ERROR: Communicator::allGatherv called with complex data\n";
} 

template <>
void
Communicator::reduce(int num, complex<double> *data, int root, MPI_Op mpi_op)
{
  cerr << "ERROR: Communicator::reduce called with complex data\n";
}

template <>
void
Communicator::broadcast(int num, complex<double> *data, int root)
{
cerr << "ERROR: Communicator::broadcast called with complex data\n";
}

#endif

