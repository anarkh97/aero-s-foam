#ifndef _CRSOLVER_H_
#define _CRSOLVER_H_

#include <Solvers.d/Solver.h>

template <class Scalar, class AnyVector, class AnyOperator>
class GenCRSolver : public GenSolver<Scalar> {
protected:
   AnyOperator *A;
   double tolerance;
   int maxiter;
   double solveTime;
 public:
   GenCRSolver(int _maxit, double _tol, AnyOperator *AA )
     { maxiter = _maxit; tolerance = _tol; A = AA; }
   void solve(AnyVector &, AnyVector &);
   int neqs() { return A->neqs(); }
   //int dim() { return A->dim(); }
   double getSolutionTime() { return solveTime; }
   long size() { return 0; }
};

#ifdef _TEMPLATE_FIX_
#include <Solvers.d/CRSolver.C>
#endif

#endif
