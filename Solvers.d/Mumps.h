#ifndef MUMPS_H_
#define MUMPS_H_

// Mumps include files
#ifdef USE_MUMPS
// include MUMPS librairies here or in header file
#include "dmumps_c.h" // double precision mumps header
#include "zmumps_c.h" // complex double precision
#endif

#include <Solvers.d/Solver.h>
#include <Math.d/SparseMatrix.h>

class EqNumberer;
class ConstrainedDSA;
template <class Scalar> class GenFullSquareMatrix;
typedef GenFullSquareMatrix<double> FullSquareMatrix;
template <class Scalar> class GenVector;
typedef GenVector<double> Vector;
template <class Scalar> class GenVectorSet;
typedef GenVectorSet<double> VectorSet;
class Connectivity;
class SolverCntl;

class FSCommunicator;

#ifdef USE_MUMPS
template<class Scalar>
class MumpsId {
};

template <>
class MumpsId<double> {
 public:
   DMUMPS_STRUC_C id;
   typedef double MumpsType;
};

template <>
class MumpsId<complex<double> > {
 public:
   ZMUMPS_STRUC_C id;
   typedef mumps_double_complex MumpsType;
};
#endif

template<class Scalar>
class GenMumpsSolver : public GenSolver<Scalar>, public GenSparseMatrix<Scalar>, public SparseData 
{
   int 	neq;        // number of equations = id.n for Mumps
   Scalar *unonz;   // matrix of elements = id.a for mumps
   int nNonZero;    // number of non zero entries = id.nz for Mumps	
   int nrbm;        // number of zero pivots detected

#ifdef USE_MUMPS
   MumpsId<Scalar> mumpsId;
#endif

   FSCommunicator *mpicomm;
   bool host;

 public:
   GenMumpsSolver(Connectivity *nToN, EqNumberer *dsa, int *map=0, FSCommunicator *_mpicomm = 0);
   GenMumpsSolver(Connectivity *nToN, DofSetArray *dsa, ConstrainedDSA *c_dsa, FSCommunicator *_mpicomm = 0);

   virtual ~GenMumpsSolver();

   void add(FullSquareMatrix &, int *dofs);
   void addImaginary(FullSquareMatrix &, int *dofs);
   void add(GenFullM<Scalar> &, int, int);
   void add(GenAssembledFullM<Scalar> &, int *);
   void addDiscreteMass(int dof, Scalar);
   void add(int dofi, int dofj, Scalar d); //HB: add upper part only 
   void addone(Scalar d, int dofi, int dofj) { add(dofi, dofj, d); }

   void unify(FSCommunicator *);
   void factor();

   void solve(Scalar *rhs, Scalar *solution);
   void reSolve(Scalar *rhs) { reSolve(1, rhs); }
   void reSolve(int nRHS, Scalar *rhs);
   void reSolve(int nRHS, Scalar **rhs);
   void reSolve(int nRHS, GenVector<Scalar> *rhs);
   void getNullSpace(Scalar *rbm);

   int dim() { return neq; }
   int neqs() { return neq; }
   long size();
   int numRBM() { return nrbm; }

   void print();
   Scalar  diag(int dof) const;
   Scalar &diag(int dof);

   void zeroAll();
   
 private:
   void init();
#ifdef USE_MUMPS
   void copyToMumpsLHS(mumps_double_complex *&m, DComplex *&d, int len);
   void copyToMumpsLHS(double *&m, double *&d, int len);
   void copyToMumpsRHS(mumps_double_complex *&m, DComplex *d, int len);
   void copyToMumpsRHS(double *&m, double *d, int len);
   void copyFromMumpsRHS(DComplex *d, mumps_double_complex *m, int len);
   void copyFromMumpsRHS(double *d, double *m, int len);
   void copyToMumpsRHS(mumps_double_complex *&m, DComplex **d, int len, int nRHS);
   void copyToMumpsRHS(double *&m, double **d, int len, int nRHS);
   void copyFromMumpsRHS(DComplex **d, mumps_double_complex *m, int len, int nRHS);
   void copyFromMumpsRHS(double **d, double *m, int len, int nRHS);
   void copyToMumpsRHS(mumps_double_complex *&m, GenVector<DComplex> *d, int len, int nRHS);
   void copyToMumpsRHS(double *&m, GenVector<double> *d, int len, int nRHS);
   void copyFromMumpsRHS(GenVector<DComplex> *d, mumps_double_complex *m, int len, int nRHS);
   void copyFromMumpsRHS(GenVector<double> *d, double *m, int len, int nRHS);
#endif
};

typedef GenMumpsSolver<double> MumpsSolver;
typedef GenMumpsSolver<DComplex> ComplexMumpsSolver;

#ifdef _TEMPLATE_FIX_
#include <Solvers.d/Mumps.C>
#endif

#endif
