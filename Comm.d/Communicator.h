#ifndef _COMMUNICATOR_H_
#define _COMMUNICATOR_H_

#include <cstdio>
#include <Utils.d/resize_array.h>
#include <Utils.d/MyComplex.h>
class Connectivity;

#ifdef USE_MPI
#ifndef MPI_NO_CPPBIND
#define MPI_NO_CPPBIND
#endif
#include <mpi.h>
#else
typedef int MPI_Comm;
typedef int MPI_Request;
typedef int MPI_Status;
typedef int MPI_Op;
typedef int MPI_Datatype;
#define MPI_COMM_WORLD 0
#define MPI_SUM 0
#endif

template <class T>
class CommTrace {
  public:
    static MPI_Datatype MPIType;
};

template <>
class CommTrace<bool> {
  private: 
    static MPI_Datatype MPIType;
};

struct RecInfo {
  int cpu, len;
};

// Utility class to contain a message link information.
struct CPair {
  int locFrom, locTo;
  int glFrom, glTo;

};

int CPairCompare(void *a, void *b);

class Communicator 
{
    MPI_Comm comm;
    ResizeArray<MPI_Request> pendReq;
    ResizeArray<MPI_Status> reqStatus;
    int nPendReq;
    int glNumCPU;
  public:
    Communicator(MPI_Comm, FILE * = stderr);
    Communicator(int _ncpu);
    Communicator(const Communicator &);
    template <class Type>
       Type globalSum(Type);
    bool globalMax(bool);
    template <class Type>
       Type globalMax(Type);
    template <class Type>
       Type globalMin(Type);
    template <class Type>
       void globalSum(int, Type*);
    template <class Type>
       void globalMax(int, Type*);
    template <class Type>
       void globalMin(int, Type*);
    template <class Type>
       void sendTo(int cpu, int tag, Type *data, int len);
    template <class Type>
       RecInfo recFrom(int tag, Type *data, int len);
    template <class Type>
       void allGather(Type *send_data, int send_count, Type *recv_data,
                      int recv_count);
    template <class Type>
       void allGatherv(Type *send_data, int send_count, Type *recv_data,
                      int recv_counts[], int displacements[]);
    template <class Type>
       void allGather(Type *recv_data, int recv_count);
    template <class Type>
       void allGatherv(Type *recv_data, int recv_counts[], int displacements[]);
    template <class Type>
       void gather(Type *send_data, int send_count, Type *recv_data,
                      int recv_count, int root);
    template <class Type>
        void gatherv(Type *send_data, int send_count, Type *recv_data,
		      int recv_counts[], int displacements[], int root);
    
    template <class Type>
      void reduce(int num, Type *data, int root = 0, MPI_Op = MPI_SUM);
    template <class Type>
      void broadcast(int num, Type* data, int root = 0);
    void sync();
    int myID();
    int numCPUs();
    MPI_Comm* getCommunicator() { return &comm; }
    void split(int, int, Communicator**);
    int remoteSize();
    void waitForAllReq(); 
};

class SysCom : public Communicator 
{
    MPI_Comm *salinasCommunicator;
  public:
    SysCom(int &argc, char **&argv);
    SysCom();
    ~SysCom();
    SysCom(MPI_Comm *communicator);
    MPI_Comm* getCommunicator() { return salinasCommunicator; }
};

// The next routine provides tag from 100 to 200 cyclicaly
int uniqueTag();

extern SysCom *syscom;
extern Communicator *structCom;

#ifdef _TEMPLATE_FIX_
#include <Comm.d/Communicator.C>
#endif

#endif
