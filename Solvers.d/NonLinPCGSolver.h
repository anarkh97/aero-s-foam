#ifndef NONLIN_PCGSOLVER_H_
#define NONLIN_PCGSOLVER_H_

class Connectivity;
class Rbm;
template <class Scalar> class GenFullSquareMatrix;
typedef GenFullSquareMatrix<double> FullSquareMatrix;

#include <Solvers.d/PCGSolver.h>

template<class Scalar, class AnyVector, class AnyOperator>
class GenNonLinPCGSolver : public GenPCGSolver<Scalar, AnyVector, AnyOperator> 
{
    int numele;
    Connectivity* allDofs;
  public:
    void reBuild(FullSquareMatrix *kel, int iter=0, int step = 1);
    void reBuild(FullSquareMatrix *kel, FullSquareMatrix *mel, double delta);

    GenNonLinPCGSolver(AnyOperator* K, int precno, int maxiter, double tolerance,
                       int numele, Connectivity* allDOFs, int kryflg,
                       int initflg, int reorthoflg, int maxVecStorage, 
                       Rbm *rbm = 0);
    virtual ~GenNonLinPCGSolver() { /* nothing to delete */ };
};

#ifdef _TEMPLATE_FIX_
#include <Solvers.d/NonLinPCGSolver.C>
#endif

#endif
