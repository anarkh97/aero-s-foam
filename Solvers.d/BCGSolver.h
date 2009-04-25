#ifndef _BCGSOLVER_H_
#define _BCGSOLVER_H_

#include <Solvers.d/Solver.h>

template <class Scalar, class AnyVector, class AnyOperator>
class GenBCGSolver : public GenSolver<Scalar> {
 protected:
   AnyOperator *A;
   double tolerance;
   int maxiter;
   double solveTime;
 public:
   GenBCGSolver(int _maxit, double _tol, AnyOperator* AA)
     { maxiter = _maxit; tolerance = _tol; A = AA; }
   ~GenBCGSolver() {};
   int neqs() { return A->neqs(); }
//   int dim() { return A->dim(); }
   void solve(AnyVector &, AnyVector &);
   double getSolutionTime() { return solveTime; }
   long size() { return 0; }
};

#ifdef _TEMPLATE_FIX_
#include <Solvers.d/BCGSolver.C>
#endif

#endif
