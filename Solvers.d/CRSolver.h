#ifndef _CRSOLVER_H_
#define _CRSOLVER_H_

#include <Solvers.d/Solver.h>

template <class Scalar, class AnyVector, class AnyOperator, class AnyPreconditioner>
class GenCRSolver : public GenSolver<Scalar> {
protected:
   AnyOperator *A;
   AnyPreconditioner *P;
   double tolerance;
   int maxiter;
   double solveTime;
 public:
   GenCRSolver(int _maxit, double _tol, AnyOperator* _A, AnyPreconditioner* __P = 0)
     { maxiter = _maxit; tolerance = _tol; A = _A; P = __P; solveTime = 0; }
   void solve(AnyVector&, AnyVector&);
   void reSolve(AnyVector &rhs) { AnyVector rhs_copy(rhs); solve(rhs_copy, rhs); }
   int neqs() { return A->neqs(); }
   double getSolutionTime() { return solveTime; }
   long size() { return 0; }
};

#ifdef _TEMPLATE_FIX_
#include <Solvers.d/CRSolver.C>
#endif

#endif
