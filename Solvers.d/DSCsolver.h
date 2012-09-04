#ifndef _DSC_SOLVER_H_
#define _DSC_SOLVER_H_

template <class Scalar> class GenVector;
typedef GenVector<double> Vector;
class Connectivity;
template <class Scalar> class GenFullM;
typedef GenFullM<double> FullM;
#include <Solvers.d/Solver.h>
#include <Math.d/SparseMatrix.h>
#ifdef DISTRIBUTED
#include <mpi.h>
#endif

class DSCsolver :
      public SparseData, public SparseMatrix, public Solver {

   double *unonz;   // storage of A

   int  *adj; 
   int *xadj;

   int scheme_number;
   int numNodes;

   double solveTime;
   double tol;

#ifdef DISTRIBUTED
   MPI_Comm dscComm;
#endif
   int color;
   int maxNum;

 public:

   DSCsolver(Connectivity *cn, EqNumberer *eqNums, int sch_number);
   virtual ~DSCsolver();

   void    add(FullSquareMatrix &knd, int *dofs) {};
   void    add(FullM &knd, int fRow, int fCol);
   void    zeroAll();

   void    factor();

   void    reSolve(double *rhs);
   void    reSolve(Vector &rhs);

   void    reSolve(int nRHS, double **RHS);
   void    reSolve(int nRHS, Vector * RHS);
   void    unify(FSCommunicator *communicator);
   void    print();

   int     dim()            { return numUncon;  }
   int     neqs()           { return numUncon;  }
   double getSolutionTime() { return solveTime; }

   // Functions that need to be written
   int numRBM() { return -1; }

   // Functions not needed
   double    diag(int) const { return 1.0; }
   double &diag(int) { throw "Crazy programmers\n"; }
   long size()    { return 0; }
};

#endif
