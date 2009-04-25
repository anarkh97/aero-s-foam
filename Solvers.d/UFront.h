#ifndef _UFRONT_H_
#define _UFRONT_H_

#include <Solvers.d/Solver.h>
#include <Utils.d/Memory.h>

#define UNROLL 6

template <class Scalar> class GenSymFullMatrix;
typedef GenSymFullMatrix<double> SymFullMatrix;
class ConstrainedDSA;
class Rbm;

class UFront : public Solver {
 protected:
  int ndof ; 		// Total number of degrees of freedom
  int maxfrsize ;	// Maximum Size of the front
  int *loc ; 		// Location of a DOF in the front. -1 if none
  int *invloc ; 	// Dof in front -1 if none
  double *data ; 	// Complete frontal data
  double **k ; 		// Pointers to row values
  double *rhs ; 	// Pointers to right hand sides
  double (*row)[UNROLL],diagBlock[UNROLL][UNROLL],factRHS[UNROLL] ;
  double (*column)[UNROLL];

  double solveTime;
  long memUsed;
  
  int rowFrontPos[UNROLL];
  int rowDOF[UNROLL];
  int currentFrontSize ;
  
  int numrbm; // number of rigid body modes
  Rbm *rbm;   // pointer to rigid body modes

  int nSaved;
  int nRows;
  int (*order)[UNROLL];
  int *rowlen;
  int (*savedFrontPos)[UNROLL];
  double (*savedDiagBlock)[UNROLL][UNROLL];
  double (*savedRHS)[UNROLL];
  double *((*savedRows)[UNROLL]) ;
  ConstrainedDSA *dsa;

  double precision; // precision for detection of RBMs

  UFront() { dsa = 0; memUsed = -memoryUsed(); }

  int locate(int i) ;
  virtual void collect(int dof);
  virtual void subFactor();
  virtual void rankNUpdate();
  virtual void saveEliminatedRows();

 public:
  UFront(int _ndof, int _frsize, double pres, Rbm *rbm = 0, int nrhs = 1);
  UFront(ConstrainedDSA* c_dsa, int _frsize, double pres, Rbm *rbm = 0,int nrhs = 1);
  ~UFront();

  virtual void finishUpdate();
  void mark_fixed(int dof) ;
  void addkel(SymFullMatrix &, int *, double*q );
  void addkel(FullSquareMatrix &, int *, double*q = 0);
  void addload(int dof, double v);
  void addmload(int dof, double *v);
  void addToDiag(int dof, double v);
  void elim(int dof) ;
  void backsub(double *q) ;
  void setDSA(ConstrainedDSA *_dsa) { dsa = _dsa; }

  int numRBM() { return numrbm; }
  void getRBMs(double *rigidBodyModes);
  void getRBMs(Vector *rigidBodyModes);
  void getRBMs(VectorSet& rigidBodyModes);

  void reSolve(double* rhs);
  void   solve(double* rhs, double *solution);
  void   solve(Vector &rhs, Vector &solution);

  // void print();
  int  neqs() { return ndof;}
  double getSolutionTime() { return solveTime; }

  long size() { return memUsed; }

};

#endif
