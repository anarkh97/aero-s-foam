#if !defined(__OPTSOL_DIST_HPP__)
#define __OPTSOL_DIST_HPP__

#if defined(STRUCTOPT) && defined(DISTRIBUTED)

#include <mpi.h>
#include <memory>
#include <Structopt.d/Optsol.h>

//------------------------------------------------------------------------------
class Optalg_dist: virtual public Optalg
{
public:
  static const int root;
  typedef enum {Optalg_begin, Optalg_end, Optalg_func, Optalg_grad} Optalg_action;

protected:
  static void bcast_action(MPI_Comm comm, Optalg_action& a);
  static void bcast_int(MPI_Comm comm, int& i);
  template<class T>
  static void bcast_array(MPI_Comm, T* ia, size_t size);
  
  static void bcast_func_params(MPI_Comm, int& iter,
				double* vars, size_t length);

  static void bcast_grad_params(MPI_Comm,
				int* active, int& actsize,
				int& funcGrad,
				double* vars, size_t varsize);
};

//------------------------------------------------------------------------------
class Optalg_dist_masterWrapper: public Optalg_dist
{
protected:
  MPI_Comm comm;
  Optalg_dist_masterWrapper(MPI_Comm c) : comm(c) {}

public:
  virtual ~Optalg_dist_masterWrapper() {}

  void solve(Optsol *_optsol);
  void func(int iter);
  void grad(int* active=0, int funcGrad=0);
};

//------------------------------------------------------------------------------
class Optalg_dist_slaveFacade: public Optalg_dist
{
private:
  MPI_Comm comm;

public:
  Optalg_dist_slaveFacade(MPI_Comm c) : comm(c) {}
  virtual ~Optalg_dist_slaveFacade() {}
  void solve(Optsol *_optsol);
  void func(int iter);
  void grad(int* active=0);

  void setOutput(FILE* f)        {}
  void buildalg(nlpdata & param) {}
  void setDefault()              {}
  void print()                   {}
  void printres()                {}
};


//------------------------------------------------------------------------------
class OptalgNlpql_dist : public OptalgNlpql, public Optalg_dist_masterWrapper
{
public:
  OptalgNlpql_dist(int subtype, MPI_Comm c) : OptalgNlpql(subtype), Optalg_dist_masterWrapper(c) {}
  virtual ~OptalgNlpql_dist() {}

  void solve(Optsol *_optsol);
  void func(int iter) { Optalg_dist_masterWrapper::func(iter); OptalgNlpql::func(iter); }
  void grad(int* active=0) { Optalg_dist_masterWrapper::grad(active); OptalgNlpql::grad(active); }

  void buildalg( nlpdata & param) { OptalgNlpql::buildalg(param); }
  void setDefault() { OptalgNlpql::setDefault(); }
  void print() { OptalgNlpql::print(); }
  void printres() { OptalgNlpql::printres(); }
};

//------------------------------------------------------------------------------
class OptalgNlpmma_dist: public OptalgNlpmma, public Optalg_dist_masterWrapper
{
public:
  OptalgNlpmma_dist(int subtype, MPI_Comm c) : OptalgNlpmma(subtype), Optalg_dist_masterWrapper(c) {}
  virtual ~OptalgNlpmma_dist() {}

  void solve(Optsol *_optsol);
  void func(int iter) { Optalg_dist_masterWrapper::func(iter); OptalgNlpmma::func(iter); }
  void grad(int* active=0) { Optalg_dist_masterWrapper::grad(active); OptalgNlpmma::grad(active); }

  void buildalg(nlpdata & param) { OptalgNlpmma::buildalg(param); }
  void setDefault() { OptalgNlpmma::setDefault(); }
  void print() { OptalgNlpmma::print(); }
  void printres() { OptalgNlpmma::printres(); }
};

/* PJSA base class OptalgSnopt is not defined
//------------------------------------------------------------------------------
class OptalgSnopt_dist: public OptalgSnopt, public Optalg_dist_masterWrapper
{
public:
  explicit OptalgSnopt_dist(MPI_Comm c) : OptalgSnopt(), Optalg_dist_masterWrapper(c) {}
  virtual ~OptalgSnopt_dist() {}

  void solve(Optsol *_optsol);
  void func(int iter) { Optalg_dist_masterWrapper::func(iter); OptalgSnopt::func(iter); }
  void grad(int* active=0) 
  { Optalg_dist_masterWrapper::grad(active, 1); OptalgSnopt::grad(active); 
    // gradient is always evaluated at the same point as function 
   }

  void buildalg(nlpdata & param) { OptalgSnopt::buildalg(param); }
  void setDefault() { OptalgSnopt::setDefault(); }
  void print() { OptalgSnopt::print(); }
  void printres() { OptalgSnopt::printres(); }  
};
*/


#endif
#endif
