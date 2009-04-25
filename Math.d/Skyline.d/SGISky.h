#ifndef _SGISKYMATRIX_H_
#define _SGISKYMATRIX_H_

#include <Math.d/SparseMatrix.h>
#include <Math.d/Skyline.d/SkyMatrix.h>

template <class Scalar> class GenSolver;
typedef GenSolver<double> Solver;
template <class Scalar> class GenFullM;
typedef GenFullM<double> FullM;
class Rbm;
class CoordSet;
class GeomState;

class SGISky : public SkyData, public SparseMatrix, public Solver 
{
   double *skyA;        // Sky array 
   double *dinv;	// inverse of diagonal

   // Timing data members
   double solveTime;
   double constructTime;

 public:
   // Constructors
   SGISky(Connectivity *cn, EqNumberer *dsa, double trbm,int *bc=0);
   SGISky(Connectivity *cn, EqNumberer *dsa, ConstrainedDSA *, 
	  double trbm, Rbm *rigid=0);
   SGISky(Connectivity *cn, EqNumberer *dsa, double trbm, int *rCN, int dummy);
   SGISky(FullM *mat, double tolerance = 1.0E-4);
   ~SGISky();

   void printMemory();
   void printConstructTime();

   double getMemoryUsed() { return 8*dlp[numUncon-1]/(1024.0*1024.0); }
   long size() { return dlp[numUncon-1]; }

   void mult(const Vector &, Vector & );     // matrix-vector multiply
   double diag( int ) const;                       // returns diagonal entry
   double &diag( int );                       // returns diagonal entry

   // assembly
   void add(FullSquareMatrix &, int *dofs);
   void add(FullM &, int rowStart, int colStart);

   void Factor();                       // Factors a matrix in sky line form.
   void Factor(Rbm *rigid);             // Factors a matrix in sky line form.

   // the following function consolidates the above Factor routines!
   void factor();
 
   void solve(double *rhs, double *solution);
   void solve(Vector &rhs, Vector &solution );

   void reSolve(Vector & rhs);               // Forward & Backward solve 
   void reSolve(double * rhs); // re solve the system, overwriting
   void reSolve(int numRHS, double **RHS);
   void reSolve(int numRHS, Vector  *RHS);

   Vector* getNullSpace();                 // retrieve ZERO ENERGY MODES 
   void    getNullSpace(double *rbm);      // retrieve ZERO ENERGY MODES 

   int  numRBM() { return 0; } // retrieve the number of rigid body modes
   void setRBM(Rbm * rigidBodyModes) { rbm = rigidBodyModes; } 
   void getRBMs(double *);       // retrieve the rigid body modes
   void getRBMs(Vector* vs);       // retrieve the rigid body modes
   void getRBMs(VectorSet& vs);       // retrieve the rigid body modes

   void addDiscreteMass(int cdof, double diMass);
   void zeroAll();
   void clean_up();
   int  dim()  { return numUncon; }
   int  neqs() { return numUncon; }
   void allMult(double x); // multiply the whole skyline matrix by x

   // Timing functions
   double getSolutionTime()  { return solveTime;     }
}; 

#endif
