#ifndef _COMM_OPT_HPP_
#define _COMM_OPT_HPP_

#ifdef STRUCTOPT

//#ifdef USE_MPI
#ifndef MPI_NO_CPPBIND
#define MPI_NO_CPPBIND
#endif

//#include <mpi.h>
#include <map>
#include <Math.d/DistVector.h>
#include <Driver.d/Communicator.h>

//------------------------------------------------------------------------------
template<class T, class Scalar>
class DistVecReducer
{
protected:
  FSCommunicator* fsComm;
  Connectivity*   nodeToSub;
  int*            glSubToLocal;
  GenSubDomain<Scalar>** subDomains;

  virtual void reduceOp(int, T**) = 0;
public:
  DistVecReducer(FSCommunicator& f, Connectivity& n2Sub, int* glS2L, GenSubDomain<Scalar>** s) : 
    fsComm(f), nodeToSub(n2Sub), glSubToLocal(glS2L), subDomains(s) {}
  virtual ~DistVecReducer() {}
  void reduce(DistVec<T>&);
};

//------------------------------------------------------------------------------
template<class Scalar>
class MaxDistVecReducer : public DistVecReducer<double, Scalar>
{
protected:
  virtual void reduceOp(int, double**);
public:
  MaxDistVecReducer(FSCommunicator& f, Connectivity& n2Sub, int* glS2L, GenSubDomain<Scalar>** s) : 
    DistVecReducer<double, Scalar>(f, n2Sub, glS2L, s) {}
};

//------------------------------------------------------------------------------
template<class T, class Scalar>
class SumDistVecReducer : public DistVecReducer<T, Scalar>
{
protected:
  virtual void reduceOp(int, T**);
public:
  SumDistVecReducer(FSCommunicator& f, Connectivity& n2Sub, int* glS2L, GenSubDomain<Scalar>** s) : 
    DistVecReducer<T, Scalar>(f, n2Sub, glS2L, s) {}
};

//------------------------------------------------------------------------------
template<class Scalar>
class PowSumDistVecReducer : public DistVecReducer<double, Scalar>
{
private:
  double pw;
protected:
  virtual void reduceOp(int, double**);
public:
  PowSumDistVecReducer(FSCommunicator& f, Connectivity& n2Sub, int* glS2L, GenSubDomain<Scalar>** s,
		       double p) : 
    DistVecReducer<double, Scalar>(f, n2Sub, glS2L, s), pw(p) {}
};

#ifdef _TEMPLATE_FIX_
#include <Structopt.d/Comm_opt.d/Comm_opt.C>
#endif

//#endif
#endif
#endif
