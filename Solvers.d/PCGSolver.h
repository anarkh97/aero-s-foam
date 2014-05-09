#ifndef _PCGSOLVER_H_
#define _PCGSOLVER_H_

#include <Solvers.d/Solver.h>
#include <Solvers.d/BasePCG.h>
#include <Solvers.d/Preconditioner.h>
#include <Solvers.d/KProject.h>
#include <Utils.d/MyComplex.h>

class Rbm;
template <class Scalar> class GenSparseMatrix;
template <class Scalar> class GenVector;

template<class Scalar, class AnyVector, class AnyOperator>
class GenPCGSolver 
	: public GenSolver<Scalar>, 
	  public BasePCG<Scalar, AnyVector, AnyOperator, KrylovProjector<Scalar,AnyVector>, Preconditioner<AnyVector> > {

protected: 

  // ... Krylov space variables
  int kryflg;
  int initflg;
  int reorthoflg;
  int maxVecStorage;

  // ... Rigid Body Mode variables
  int numrbm;		// # of rigid body modes
  Rbm *rbm;		// pointer to rigid body modes

  //double solveTime;
  //long memUsed;

public:
  // ... Constructors
  GenPCGSolver(AnyOperator *A, int precno, int maxiter, double tolerance, int maxVecStorage = 0,
               Rbm *rbm=0);
  // destructor
  virtual ~GenPCGSolver() { /* nothing to delete */ };

  // ... Rigid Body Mode functions
  int  numRBM() { return numrbm; }
  void getRBMs(double *rigidBodyModes);
  void getRBMs(Vector *rigidBodyModes);
  void getRBMs(VectorSet &rigidBodyModes);

  // ... Linear solution functions
  void solve(Scalar* rhs, Scalar* solution); 
  void solve(AnyVector &rhs, AnyVector &solution); 

  // ... Linear solution functions that overwrite the rhs vector
  void reSolve(Scalar* rhs);  
  void reSolve(AnyVector &rhs);

  // Multiple rhs reSolve functions
  void reSolve(int nRHS, Scalar **RHS);  
  void reSolve(int nRHS, AnyVector *RHS);

  long size() { return 0; }
  int neqs(); 

  void factor() { };

};

typedef GenPCGSolver<double, GenVector<double>, GenSparseMatrix<double> > PCGSolver;

#ifdef _TEMPLATE_FIX_
  #include <Solvers.d/PCGSolver.C>
#endif

#endif

