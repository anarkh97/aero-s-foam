#ifndef _BLOCK_SKY_H_
#define _BLOCK_SKY_H_

/*****************************************************************************
 *                   Copyright (C) 1999 CMSoft                               *
 *                                                                           *
 *  These lines of code and declarations contain unpublished proprietary     *
 *  information of CMSoft. They may not be copied or duplicated in whole     *
 *  or part without prior authorization from CMSoft.                         *
 *                                                                           *
 *****************************************************************************/

#ifdef sgi
#include <ulocks.h>
#endif

class Connectivity;
class EqNumberer;
class DofSetArray;
template <class Scalar> class GenFullM;
typedef GenFullM<double> FullM;
template <class Scalar> class GenAssembledFullM;
typedef GenAssembledFullM<double> AssembledFullM;
template <class Scalar> class GenVector;
typedef GenVector<double> Vector;

#include <cstdio>
#include <Math.d/SparseMatrix.h>
#include <Solvers.d/Solver.h>
#include <Utils.d/MyComplex.h>

template<class Scalar>
class GenBlockSky : public GenSolver<Scalar>, public GenSparseMatrix<Scalar> 
{
     int nBlocks;   // number of blocks or supernodes
     int *blWeight; // multiplicity of each supernode
     int *blHeight; // height of each supernodal rectangle
     int *blTop; // begining of each supernodal rectangle
     int *firstDof;
     int *lastCol;
     int *dlp; // diagonal location pointer (non-FORTRAN style)
     int myRCN;
     int *rowColNum;
     int maxBlockSize; // maximum area of a block
     Scalar *invDiag;
     int *perm;

     double tol;

     int neq;
     int nzem;
     int *sing; // PJSA: location of singularities (used for GtGstar in FetiDPSolver)
   protected:
     Scalar *skyA;
     #ifdef sgi
       void pFactor(int iThread, int numThread, barrier_t *b, 
                    Scalar *ltmp, Scalar *invDiag, Scalar *origDiag);
     #endif
   public:
     GenBlockSky(Connectivity *nodeToNode, EqNumberer *eqnums, double tol);
     GenBlockSky(Connectivity *nodeToNode, DofSetArray *eqnums, double tol);
     GenBlockSky(Connectivity *nodeToNode, DofSetArray *dsa, double tol, 
              int *glInternalMap);
     virtual ~GenBlockSky();
     void factor();
     void unify(FSCommunicator *communicator);
     void parallelFactor();
     void solve(Scalar *rhs, Scalar *solution);
     void solve(GenVector<Scalar> &rhs, GenVector<Scalar> &solution );
     void reSolve(Scalar *f);
     void reSolve(GenVector<Scalar> &f);
     void reSolve(int numRHS, Scalar **RHS);
     // assembly
     void add(FullSquareMatrix &, int *dofs);
     void add(FullSquareMatrixC &, int *dofs);
     void addImaginary(FullSquareMatrix &kel, int *dofs);
     void add(AssembledFullM &, int *dofs);
     void add(GenAssembledFullM<complex<double> > &, int *dofs);
     void add(FullM &, int rowStart, int colStart);
     void add(int row_dof, int col_dof, Scalar s);
     void addBoeing(int, const int *, const int *, const double *, int *, Scalar multiplier);
     void addDiscreteMass(int dof, Scalar dmass);

     void print(FILE * = stderr);
     int neqs() { return neq; }
     double getSolutionTime() { return 0.0; }
     long size() { return dlp[neq-1]; }
     Scalar diag(int i) const { return skyA[dlp[i]]; }  // Jing Li's problem
     Scalar &diag(int i) { return skyA[dlp[i]]; }  // Jing Li's problem
     void zeroAll();
     void clean_up();
     int dim()    { return neq;  }
     int numRBM() { return nzem; }
     int* getSingularities() { return sing; }
   private:
     void initialize();
};

typedef GenBlockSky<double> BlockSky;

#endif
