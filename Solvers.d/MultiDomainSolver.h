#ifndef MULTI_DOMAIN_SOLVER
#define MULTI_DOMAIN_SOLVER

#include <Solvers.d/ParallelSolver.h>

template<class Scalar> class GenDistrVector;
template<class Scalar> class GenSubDomain;
class FSCommunicator;

template <class Scalar> 
class MultiDomainSolver : public GenParallelSolver<Scalar>
{
  protected:
    int neq;
    int nsub;
    GenSubDomain<Scalar> **sd;
    FSCommunicator *com;
  public:
    MultiDomainSolver() {}
    MultiDomainSolver(int _neq, int _nsub, GenSubDomain<Scalar> **_sd, FSCommunicator *_com) : neq(_neq), nsub(_nsub), sd(_sd), com(_com) {}
    void reSolve(GenDistrVector<Scalar> &);
    void solve(GenDistrVector<Scalar> &, GenDistrVector<Scalar> &);
    virtual void solve(Scalar *rhs, Scalar *solution) = 0;
    double getFNormSq(GenDistrVector<Scalar> &f);
  private:
    void multLTinv(Scalar *f_g, GenDistrVector<Scalar> &f);
    void multL(Scalar *u_g, GenDistrVector<Scalar> &u);
    void multLinv(GenDistrVector<Scalar> &u, Scalar *u_g);
    void multLT(GenDistrVector<Scalar> &f, Scalar *f_g);
};

#ifdef _TEMPLATE_FIX_
#include <Solvers.d/MultiDomainSolver.C>
#endif

#endif
