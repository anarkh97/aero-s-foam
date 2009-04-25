#ifndef _SGISPARSEMATRIX_H_
#define _SGISPARSEMATRIX_H_

template <class Scalar> class GenVector;
typedef GenVector<double> Vector;
class Connectivity;
class DofSetArray;
class ConstrainedDSA;
template <class Scalar> class GenFullSquareMatrix;
class Rbm;
template <class Scalar> class GenSolver;

#include <Math.d/SparseMatrix.h>
#include <Utils.d/MyComplex.h>

template<class Scalar>
class GenSGISparseMatrix :
        public SparseData, public GenSparseMatrix<Scalar>, public GenSolver<Scalar> {

protected:
   Scalar *unonz;
   Scalar *newrhs;
   double ops;
   int nonz;
   int  token;
   Rbm *rbm;
   int  numrbm;
   double solveTime;
   double constructTime;

 public:
   GenSGISparseMatrix(Connectivity *, DofSetArray *, ConstrainedDSA *, Rbm *rbm=0,
                      int token=0);
   virtual ~GenSGISparseMatrix();

   Scalar  diag(int dof) const;
   Scalar &diag(int dof);
   void    add(FullSquareMatrix &, int *dofs);
   void    zeroAll();
   void    clean_up();
   void    print(const char *msg = "");
   void    factor();
   void    solve(Scalar *rhs, Scalar *solution);
   void    solve(GenVector<Scalar> &rhs, GenVector<Scalar> &solution);
   void    reSolve(Scalar *rhs);
   void    reSolve(GenVector<Scalar> &rhs);
   void    reSolve(int nRHS, Scalar **RHS);
   void    reSolve(int nRHS, GenVector<Scalar> *RHS);
   int     numRBM()  { return numrbm; }
   void    getRBMs(double *);    
   void    getRBMs(Vector *);   
   void    getRBMs(VectorSet &);
   void    addDiscreteMass(int dof, Scalar mass);
   double  getMemoryUsed();
   long size();
   int     dim()              { return numUncon;      }
   int     neqs()             { return numUncon;      }
   double  getSolutionTime()  { return solveTime;     }
   double  getConstructTime() { return constructTime; }
};

typedef GenSGISparseMatrix<double> SGISparseMatrix;
typedef GenSGISparseMatrix<DComplex> ComplexSGISparseMatrix;

#ifdef _TEMPLATE_FIX_
#include <Math.d/SGISparseMatrix.C>
#endif

#endif
