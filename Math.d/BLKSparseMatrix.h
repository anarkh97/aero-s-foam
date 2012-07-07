#ifndef BLOCK_SPARSE_MATRIX_H_
#define BLOCK_SPARSE_MATRIX_H_

template <class Scalar> class GenVector;
typedef GenVector<double> Vector;
template <class Scalar> class GenVectorSet;
typedef GenVectorSet<double> VectorSet;
class Rbm;
template <class Scalar> class GenSparseMatrix;
typedef GenSparseMatrix<double> SparseMatrix;
template <class Scalar> class GenSolver;
typedef GenSolver<double> Solver;
class DofSetArray;
class DofSetArray;
class Connectivity;
template <class Scalar> class GenFullM;
typedef GenFullM<double> FullM;

#include <Math.d/SparseMatrix.h>
#include <Utils.d/MyComplex.h> 
#include <Solvers.d/Solver.h>

template<class Scalar>
class GenBLKSparseMatrix :
        public SparseData, public GenSparseMatrix<Scalar>, public GenSolver<Scalar> {

protected:

   Scalar *unonz;

   int nsuper;
   int nsub; 
   int nnzl; 
   int tmpsiz;

   int *colcnt;
   int *snode;
   int *xsuper;
   int *invsuper;

   int *adj; 
   int *xadj;

   int *iwork;
   int iwsiz;

   int *xlindx;
   int *lindx;

   int  *xlnz;
   Scalar *lnz;

   int *perm;
   int *invp;

   double solveTime;
   double tol;
   
   Rbm *rbm;   // pointer to rigid body modes
   int numrbm; // number of rigid body modes
   int ngrbm;  // number of geometric rbms
   int defblk; // size of last deficient block
   int lbdef; 
   int *iprow;
   int *ipcol;
   int *def;

   int spRenum; // 0 = esmond, 1 = metis
   bool myRbm;

 public:

   // Stiffness matrix
   GenBLKSparseMatrix(Connectivity *, DofSetArray *, DofSetArray *,
                      double tolerance, int spRenum = 0, Rbm *rbm=0);
   // Kii
   GenBLKSparseMatrix(Connectivity *, DofSetArray *, int *dofmap,
                      double tolerance, int spRenum = 0, Rbm *rbm=0);
   // GtG, Kcc
   GenBLKSparseMatrix(Connectivity *, EqNumberer *, double tolerance, int spRenum = 0, int ngrbm = 0);

   virtual ~GenBLKSparseMatrix();

   void allocateMemory();

   Scalar  diag(int dof) const;
   Scalar  &diag(int dof);
   void    unify(FSCommunicator *communicator);
   void    addBoeing(int nlines, const int *Kai, const int *Kaj,
                     const double *nz, int *map, Scalar multiplier);
   void    add(FullSquareMatrix &, int *dofs);
   void    add(FullSquareMatrixC &, int *dofs);
   void    add(FullM &knd, int fRow, int fCol);
   void    add(GenAssembledFullM<Scalar> &kel, int *dofs);
   void    add(int row_dof, int col_dof, Scalar s) { addone(s, row_dof, col_dof); }
   void    addone(Scalar d, int dofi, int dofj);
   void    add(Scalar *_lnz);
   Scalar  getone(int row, int col);
   void    zeroAll();
   void    clean_up();
   int     dim()  { return numUncon; }
   int     neqs() { return numUncon; }
   void    printAll();
   void    print();
   Scalar* getData() { return lnz; }

   void    factor();

   void    solve(Scalar *rhs, Scalar *solution);
   void    solve(GenVector<Scalar> &rhs, GenVector<Scalar> &solution );

   void    reSolve(Scalar *rhs);
   void    reSolve(GenVector<Scalar> &rhs);

   void    reSolve(int nRHS, Scalar **RHS);
   void    reSolve(int nRHS, GenVector<Scalar> * RHS);

   void    reSolve(GenFullM<Scalar> *);

   double getMemoryUsed();
   long size() { return (numUncon) ? xlnz[numUncon] : 0; }

   void    getRBMs(double *);
   void    getRBMs(Vector *);
   void    getRBMs(VectorSet &);
   int     numRBM();

   void    addDiscreteMass(int dof, Scalar mass);
   void    addImaginary(FullSquareMatrix &kel, int *dofs);
   double getSolutionTime() { return solveTime; }

   void mult(const Scalar *rhs, Scalar *result);
   void getNullSpace(Scalar *ns);

 private:
   void init();
   void computeRBMs();
};

template<class Scalar>
class WrapSparseMat : public GenBLKSparseMatrix<Scalar>
{
  public:
    struct CtorData {
      Connectivity *cn;
      DofSetArray *dsa, *cdsa;
      double trbm;
      Rbm *rbm;
      int spRenum;
      CtorData(Connectivity *c, DofSetArray *d, DofSetArray *dc, double t, int sr, Rbm *r) {
        cn = c;
        dsa = d;
        cdsa = dc;
        trbm = t;
        spRenum = sr;
        rbm = r;
      }
    };

    WrapSparseMat(CtorData &ctd) : GenBLKSparseMatrix<Scalar>(ctd.cn, ctd.dsa, ctd.cdsa,
        ctd.trbm, ctd.spRenum, ctd.rbm) {}
};

typedef GenBLKSparseMatrix<double> BLKSparseMatrix;
typedef GenBLKSparseMatrix<DComplex> BLKSparseMatrixC;

#ifdef _TEMPLATE_FIX_
#include <Math.d/BLKSparseMatrix.C>
#endif

#endif
