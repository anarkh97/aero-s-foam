#ifndef _FETIOP_H_
#define _FETIOP_H_

#include <Threads.d/Paral.h>

template <class Scalar> class GenFetiOpControler;
template <class Scalar> class GenFetiSolver;
class IntFullM;
template <class Scalar> class GenSubDomain;
template <class Scalar> class GenSparseMatrix;

template<class Scalar>
class GenFetiOp : public TaskDescr 
{
   protected:
    GenSubDomain<Scalar> *sd;
    GenSolver<Scalar> *solver;
    GenSolver<Scalar> *K;
    GenSparseMatrix<Scalar> *KasSparse;
    Rbm *rbm; // For dynamics only and nonlinear
    int numNeighb;
    int halfOffset;
    int *alphaOffset;
    int *betaOffset;
    int isFeti2, isDynamic;

    int numRBM;            // total number of RBMs
    double *locRBMs;       // local RBMs for the whole domain
    double *locInterfRBMs; // local RBMs on the interface in local linear form
                           // The modes are consecutive in the array

    int *neighbNumRBMs;
    double **neighbRBMs;

    int crnDofSize;        // number of additional lagrange multipliers
    IntFullM *BClocal;

    GenFetiOpControler<Scalar> *control;

    int QGisLocal;  // Wether or not QG is local

    Scalar *interfBuff;
    FSCommPattern<Scalar> *vPat;

   public:
    GenFetiOp() { init(); }
    GenFetiOp(GenSubDomain<Scalar> *, GenFetiOpControler<Scalar> *, 
              int, int, FSCommPattern<Scalar> *, Rbm * =0);
    virtual ~GenFetiOp();
    void setSysMatrix(GenSolver<Scalar> *k, GenSparseMatrix<Scalar> *ks) { K = k; KasSparse = ks; }
    void run() override;
    void runFor(int) override { throw "Illegal operation called on GenFetiOp"; }
    void clean_up();

    void localSolve();
    void sendInterfRBM(FSCommPattern<Scalar> *rbmPat);
    void sendNumNeighbRBM(FSCommPattern<int> *sPat);
    void getNumNeighbRBM(FSCommPattern<int> *sPat);
    void getNeighbQGs(FSCommPattern<Scalar> *rbmPat);
    void getGtMult();
    void getGtQMult();
    void getGtFMult();
    void subAlphaG();
    void subAlphaG1();
    void subAlphaG2();
    void subAlphaGQ();
    void subNuC();
    void initializeCRNs(FSCommPattern<int> *sPat);
    void assembleGtCs();
    void getNumNeighbCRNs(FSCommPattern<int> *sPat);
    void assembleGtQGs();
    void makeCoarseSet();
    void computeFiBC();
    void getNeighbFC();
    void assembleGtFCs();
    void assembleCtFCs();
    void computeNeighborFGs();
    void getCtFMult();
    void getCtMult();
    void subNuFC();
    void subAlphaFG();
    void setHalfOffset(int a) { halfOffset = a; }
    void reSetAlphaOffsets(int *v);
    void setAlphaOffsets(int *v);
    void setBetaOffsets(int *v);

    int  getNumRBM()     { return numRBM;      }
    int  getcrnDofSize() { return crnDofSize;  }
    int  numEdge()       { return numNeighb;   } 

    void setglobalSum(void *);
    double res;
    friend class GenFetiSolver<Scalar>;

  private:
    void init();
};

typedef GenFetiOp<double> FetiOp;

#ifdef _TEMPLATE_FIX_
  #include <Feti.d/FetiOp.C>
#endif

#endif


