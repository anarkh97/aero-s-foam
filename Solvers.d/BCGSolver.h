#ifndef _BCGSOLVER_H_
#define _BCGSOLVER_H_

#include <Solvers.d/Solver.h>

template <class Scalar, class AnyVector, class AnyOperator, class AnyPreconditioner>
class GenBCGSolver : public GenSolver<Scalar> {
 protected:
   AnyOperator *A;
   AnyPreconditioner *P;
   double tolerance;
   int maxiter;
   double solveTime;
 public:
   GenBCGSolver(int _maxit, double _tol, AnyOperator *_A, AnyPreconditioner *_P = 0)
     { maxiter = _maxit; tolerance = _tol; A = _A; P = _P; }
   ~GenBCGSolver() {};
   int neqs() { return A->neqs(); }
   void solve(AnyVector &, AnyVector &);
   double getSolutionTime() { return solveTime; }
   long size() { return 0; }
};

#ifdef _TEMPLATE_FIX_
#include <Solvers.d/BCGSolver.C>
#endif

#endif
