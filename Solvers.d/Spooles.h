#ifndef SPOOLES_H_
#define SPOOLES_H_

// Spooles include files
#ifdef USE_SPOOLES
extern "C" {
  #include <misc.h>
  #include <FrontMtx.h>
  #include <SymbFac.h>
}
#endif

#include <Solvers.d/Solver.h>
#include <Math.d/SparseMatrix.h>
#include <Utils.d/MyComplex.h>


class ConstrainedDSA;
template <class Scalar> class GenFullSquareMatrix;
typedef GenFullSquareMatrix<double> FullSquareMatrix;
typedef GenFullSquareMatrix<DComplex> FullSquareMatrixC;

template <class Scalar> class GenVector;
typedef GenVector<double> Vector;
template <class Scalar> class GenVectorSet;
typedef GenVectorSet<double> VectorSet;
class Connectivity;

class FSCommunicator;

template<class Scalar>
class GenSpoolesSolver : public GenSolver<Scalar>, public GenSparseMatrix<Scalar>,
                         public SparseData 
{
   int neq;               // number of equations
   int *constrNum;        // constrained equation numbers
   Scalar *unonz;
   double cpus[22];
   int stats[7];
   int pivotingflag;
   int numThreads;
   int nNonZero;
   int isScaled;        // whether to scale the matrix or not
   Scalar *scale;       // vector to store the matrix scaling
   int msglvl;
   FILE *msgfile;
   long _size;

#ifdef USE_SPOOLES
   InpMtx   *inpMtx;	 // input Matrix
   FrontMtx *frontMtx;    // front Matrix
   IV       *newToOldIV;  // new to old permutation table
   IV       *oldToNewIV;  // old to new permutation table
   SubMtxManager *mtxManager;
   IV *ownersIV;
   DenseMtx *mtxB, *mtxX;
   IVL *symbfacIVL; // PJSA
   ETree *frontETree; // PJSA
   DV *cumopsDV; // PJSA
   Graph *graph; // PJSA
#endif

 public:
   //GenSpoolesSolver(Connectivity *nToN, DofSetArray *dsa, int *map=0);
   GenSpoolesSolver(Connectivity *nToN, EqNumberer *dsa, int *map=0);
   GenSpoolesSolver(Connectivity *nToN, DofSetArray *_dsa, ConstrainedDSA *c_dsa);

   virtual void clean_up() {
     cleanUp();
     if(unonz) { delete [] unonz; unonz = 0; }
     if(scale) { delete [] scale; scale = 0; }
   }

   virtual ~GenSpoolesSolver();

   void add(FullSquareMatrix &, int *dofs);
   void addImaginary(FullSquareMatrix &, int *dofs);
   void add(FullSquareMatrixC&, int *dofs); // RT addded to support PML, DGM
   void add(GenFullM<Scalar> &, int *dofs);
   void add(GenFullM<Scalar> &, int, int);
   void add(GenAssembledFullM<Scalar> &, int *);
   void addDiscreteMass(int dof, Scalar);
   void add(int dofi, int dofj, Scalar d); //HB: add upper part only 
   void addone(Scalar d, int dofi, int dofj) { add(dofi, dofj, d); }
   //void addBoeing(int nlines, const int *Kai, const int *Kaj,
   //               const double *nz, int *map, Scalar multiplier);

   void unify(FSCommunicator *);
   void parallelFactor();
   void factor();
   void allFactor(bool fctIsParal);

   void solve(Scalar *rhs);
   void reSolve(Scalar *rhs) { solve(rhs); }
   void solve(Scalar *rhs, Scalar *solution);
   
   void print();
   int dim()  { return neq;      }
   int neqs() { return neq; }
   double getMemoryUsed();
   long size();
   Scalar  diag(int dof) const;
   Scalar &diag(int dof);

   void    zeroAll();
   void    cleanUp();
   double  getSolutionTime()  { return 0.0; }
   double  getConstructTime() { return 0.0; }

   int numRBM() { return 0; } // note: spooles should not be used for singular matrices.
 
 private:
   void init();
   void applyScaling(Scalar *vector);
   void symmetricScaling();
};

template<class Scalar>
class WrapSpooles : public GenSpoolesSolver<Scalar>
{
  public:
    struct CtorData {
      Connectivity *cn;
      DofSetArray *dsa;
      ConstrainedDSA *cdsa;
      Rbm *rbm;
      CtorData(Connectivity *c, DofSetArray *d, ConstrainedDSA *dc) {
        cn = c;
        dsa = d;
        cdsa = dc;
      }
    };

    WrapSpooles(CtorData &ctd) : GenSpoolesSolver<Scalar>(ctd.cn, ctd.dsa, ctd.cdsa) {}
};

typedef GenSpoolesSolver<double> SpoolesSolver;
typedef GenSpoolesSolver<DComplex> ComplexSpoolesSolver;

#endif
