#ifndef _SKYMATRIX_H_
#define _SKYMATRIX_H_

#include <cstdio>

#include <Math.d/SparseMatrix.h>
#include <Solvers.d/Solver.h>
#include <Utils.d/MyComplex.h>
#if defined(sgi) && ! defined(_OPENMP)
#include <ulocks.h>
#endif

class Rbm;
template <class Scalar> class GenFullSquareMatrix;
typedef GenFullSquareMatrix<double> FullSquareMatrix;
typedef GenFullSquareMatrix<DComplex> FullSquareMatrixC;
class CoordSet;

class SkyData {
 protected:
   int *lacol;          // last active column (lacol)
   int *pivot;          // stores pivot information for the solver
   int *seqid;          // stores pivot information for geometric sky solver
   int *dlp;            // diagonal location pointer

                        // in some instances, rowColNum is allocated
                        // from a SkyData constructor and other times,
                        // it is allocated in the ConstrainedDSA class.
   int myRCN;           // 0=not allocated, 1=allocated memory
   int *rowColNum;      // unconstrained row/col numbers

   int neq;             // # of equations
   int numUncon;        // # unconstrained dofs
   double TOLERANCE;    // Tolerance value for factoring

   int nzem;            // # zero energy modes
   Rbm *rbm;            // stores rigid body mode information
   int isTRBM;          // whether we are TRBM or GRBM

   int myRbm;

 public:
   SkyData();
   SkyData(Connectivity *cn, EqNumberer *dsa, double trbm, int *bc);
   SkyData(Connectivity *cn, EqNumberer *dsa, double trbm);
   SkyData(Connectivity *cn, DofSetArray *c_dsa, double trbm, Rbm *rigid=0);
   SkyData(EqNumberer *_dsa, Connectivity *cn, double trbm, int *rCN);
   SkyData(int n, double tolerance = 1.0E-4);
   virtual ~SkyData();
 private:
   void initialize();
};

// Sky Matrix:

template<class Scalar>
class GenSkyMatrix : public SkyData, public GenSparseMatrix<Scalar>, public GenSolver<Scalar>
{
protected:
   Scalar *skyA;        // Sky array
   int isScaled;        // whether to scale the matrix or not.
   Scalar *scale;       // vector to store the matrix scaling
   bool wasScaled;      // HB: indicate if the scaling has already been applied before
                        //     getting in Factor/parallelFactor method
   // Timing data members
   double solveTime;
   double constructTime;
 public:
   // Constructors
   GenSkyMatrix() { this->skyZ = 0; }
   GenSkyMatrix(Connectivity *cn, EqNumberer *dsa, double trbm, int *bc);
   GenSkyMatrix(Connectivity *cn, EqNumberer *dsa, double trbm, int isScaled=0);
   GenSkyMatrix(Connectivity *cn, DofSetArray *, double trbm, Rbm *rigid=0);
   GenSkyMatrix(Connectivity *cn, EqNumberer *dsa, double trbm, int *rCN, int dummy);
   GenSkyMatrix(GenFullM<Scalar> *mat, double tolerance = 1.0E-4);
   GenSkyMatrix(Connectivity *, EqNumberer *, ConstrainedDSA *, double trbm);
   GenSkyMatrix(int, double);
   // Destructor
   virtual ~GenSkyMatrix();

   // returns memory used in megabytes
   double getMemoryUsed() { return 8*dlp[numUncon-1]/(1024.0*1024.0); }
   long size() { return (numUncon) ? dlp[numUncon - 1] : 0; }
   Scalar* getData() { return skyA; }

   void mult(const GenVector<Scalar> &, GenVector<Scalar> &);
   void mult(const Scalar *rhs, Scalar *result);
   Scalar diag(int) const;                           // returns diagonal entry
   Scalar &diag(int);

   // assembly
   void add(FullSquareMatrix &, int *dofs);
   void add(FullSquareMatrixC &, int *dofs);
   void addImaginary(FullSquareMatrix &, int *dofs);
   void add(GenAssembledFullM<Scalar> &kel, int *dofs);
   void add(GenFullM<Scalar> &, int rowStart, int colStart);
   void addone(Scalar d, int dofi, int dofj);
   void addBoeing(int, const int *, const int *, const double *, int *, Scalar multiplier);
   void add(Scalar *_skyA);
   void addPoint(Scalar, int, int);

   void Factor();
   void Factor(Rbm *rigid);
   void factor();
   // Parallel factorization. It's used for Coarse problems only
   virtual void parallelFactor();

   void solve(Scalar *rhs, Scalar *solution);
   void solve(GenVector<Scalar> &rhs, GenVector<Scalar> &solution);

   void reSolve(GenVector<Scalar> &rhs);
   void reSolve(Scalar *rhs);
   void reSolve(int numRHS, Scalar **RHS);
   void reSolve(int numRHS, GenVector<Scalar> *RHS);
   void reSolve(int numRHS, Scalar *RHS);
   void reSolve(GenFullM<Scalar> *mat);

   void forward(GenVector<Scalar> &rhs);
   void backward(GenVector<Scalar> &rhs);

   void unify(FSCommunicator *communicator);

   GenVector<Scalar>* getNullSpace();   // retrieve ZERO ENERGY MODES
   void getNullSpace(Scalar *rbm);      // retrieve ZERO ENERGY MODES

   int  numRBM() { return nzem; } // retrieve the number of rigid body modes
   void setRBM(Rbm * rigidBodyModes) { rbm = rigidBodyModes; }
   void getRBMs(double *);         // retrieve the rigid body modes
   void getRBMs(Vector* vs);       // retrieve the rigid body modes
   void getRBMs(VectorSet& vs);    // retrieve the rigid body modes

   void addDiscreteMass(int cdof, Scalar diMass);
   void add(int row_dof, int col_dof, Scalar s) { addone(s, row_dof, col_dof); }
   void zeroAll();
   int  dim()  { return numUncon; }
   int  neqs() { return numUncon; }
   void allMult(Scalar x); // multiply the whole skyline matrix by x

   // Timing functions
   double getSolutionTime()  { return solveTime; }

   void printMemory();
   void printConstructTime();
   void print1(int dof);
   void print(FILE * =stderr);
   void printMatlab(int subNumber);
   void printMatlab(char *fileName);
   void printDiagonals();

   void clean_up();

   Scalar getone(int row, int col);

   void applyScaling(Scalar *vector);
   void symmetricScaling();
   bool IsScaled() { return (isScaled)? true : false; }
   void setIsScaled(int _isScaled) { isScaled = _isScaled; }

   //HB
   double rmsBandwidth();

#if defined(sgi) && ! defined(_OPENMP)
   void pfact(int, int, barrier_t *, Scalar *);
#else
   void pfact(int, int, Scalar *);
#endif
};

template<class Scalar>
class WrapSkyMat : public GenSkyMatrix<Scalar>
{
  public:
    struct CtorData {
      Connectivity *cn;
      DofSetArray *dsa;
      double trbm;
      Rbm *rbm;
      CtorData(Connectivity *c, DofSetArray *d, double t, Rbm *r) {
        cn = c;
        dsa = d;
        trbm = t;
        rbm = r;
      }
    };

    WrapSkyMat(CtorData &ctd) : GenSkyMatrix<Scalar>(ctd.cn, ctd.dsa, ctd.trbm, ctd.rbm) {}
};


typedef GenSkyMatrix<double> SkyMatrix;
typedef GenSkyMatrix<DComplex> SkyMatrixC;

#ifdef _TEMPLATE_FIX_
#include <Math.d/Skyline.d/SkyMatrix.C>
#endif


#endif
