#ifndef __GMRESSOLVER_H__
#define __GMRESSOLVER_H__

#include <Solvers.d/Solver.h>

template <class Scalar> class GmresOrthoSet;
class FSCommunicator;

template<class Scalar, class AnyVector, class AnyOperator, class LeftPreconditioner, class RightPreconditioner>
class GmresSolver : public GenSolver<Scalar>
{
  private: 
    int maxit, maxortho, printNumber, verbose;
    double tol;
    AnyOperator *op;
    void (AnyOperator::*matvec)(AnyVector &, AnyVector &);
    LeftPreconditioner *leftprec;
    void (LeftPreconditioner::*applyLeft)(AnyVector &, AnyVector &);
    RightPreconditioner *rightprec;
    void (RightPreconditioner::*applyRight)(AnyVector &, AnyVector &);
    GmresOrthoSet<Scalar> *oSetGMRES;
    FSCommunicator *com; 
    int rank;
    int m_info[1];

  public:
    // Constructor
    GmresSolver(int maxit, double tol, AnyOperator *_op, void (AnyOperator::*_matvec)(AnyVector &, AnyVector &),
                LeftPreconditioner *_leftprec, void (LeftPreconditioner::*_applyLeft)(AnyVector &, AnyVector &), 
                RightPreconditioner *_rightprec, void (RightPreconditioner::*_applyRight)(AnyVector &, AnyVector &),
                FSCommunicator* _comm = NULL);

    // Destructor
    ~GmresSolver();

    // Linear solution function
    void solve(AnyVector &b, AnyVector &x) { x=b; reSolve(x); }
    void reSolve(AnyVector &x); 
    void reset() { oSetGMRES->reset(); }
    int info(int i) { return m_info[i]; }
    long size() { return 0; }
    int neqs() { return op->neqs(); }
    void factor() { }
};

#ifdef _TEMPLATE_FIX_
  #include <Solvers.d/GmresSolver.C>
#endif

#endif
