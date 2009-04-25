#ifndef BLOCK_SPARSE_MATRIXC_H_
#define BLOCK_SPARSE_MATRIXC_H_

class ComplexVector;
class ComplexVectorSet;
class DofSetArray;
class ConstrainedDSA;
class Connectivity;
class FullMC;
class AssembledFullMC;

#include <Utils.d/MyComplex.h>
#include <Math.d/SparseMatrix.h>
#include <Solvers.d/ComplexSolver.h>

class BLKSparseMatrixC :
        public SparseData, public SparseMatrix, public ComplexSolver {

protected:

   DComplex *unonz;

   int nsuper;
   int   nsub; 
   int   nnzl; 
   int tmpsiz;

   int *colcnt;
   int *snode;
   int *xsuper;
   int *invsuper;

   int  *adj; 
   int *xadj;

   int *iwork;
   int  iwsiz;

   int *xlindx;
   int *lindx;

   int    *xlnz;
   DComplex *lnz;

   int *perm;
   int *invp;

   double solveTime;
   double tol;

   int numrbm; // number of rigid body modes
   int defblk; // size of last deficient block
   int lbdef;  
   int *iprow; // row pivoting sequence for last block
   int *ipcol; // col pivoting sequence for last block
   int *def;   // indices of columns linearly dependent

 public:

   // Stiffness matrix
   BLKSparseMatrixC(Connectivity *, DofSetArray *, ConstrainedDSA *,
                   double tolerance);
   // Kii
   BLKSparseMatrixC(Connectivity *, DofSetArray *, int *dofmap,
                   double tolerance);
   // GtG 
   BLKSparseMatrixC(Connectivity *, EqNumberer *, double tolerance);

   void allocateMemory();

   double diag(int dof) {return 0.0;}
   DComplex diagComplex(int dof);

   void    add(FullSquareMatrix &, int *dofs);
   void    add(FullM &knd, int fRow, int fCol);

   void    add(FullSquareMatrixC &, int *dofs);
   void    add(AssembledFullMC &, int *dofs);
   void    add(FullMC &knd, int fRow, int fCol);
   void    addImaginary(FullSquareMatrix &, int *dofs);
   void    addImaginary(FullM &knd, int fRow, int fCol);

   void    zeroAll();
   int     dim()  { return numUncon; }
   int     neqs() { return numUncon; }
   void    printAll();

   void    factor();

   void    solve(DComplex *rhs, DComplex *solution);
   void    solve(ComplexVector &rhs, ComplexVector &solution );

   void    reSolve(DComplex *rhs);
   void    reSolve(ComplexVector &rhs);

   void    reSolve(int nRHS, DComplex *RHS);
   void    reSolve(int nRHS, DComplex **RHS);
   void    reSolve (FullMC *);


   double getMemoryUsed();
   long size() { return xlnz[numUncon]; }

   int     numRBM() { return numrbm; };

   double getSolutionTime() { return solveTime; }

};

#endif
