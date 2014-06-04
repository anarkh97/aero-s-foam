#ifndef _STATIC_DESCR_H_
#define _STATIC_DESCR_H_

#include <Utils.d/MyComplex.h>

class Domain;
template <class Scalar> class GenSolver;
typedef GenSolver<double> Solver;
typedef GenSolver<DComplex> ComplexSolver;
template <class Scalar> class GenSparseMatrix;
typedef GenSparseMatrix<double> SparseMatrix;
class GeomState;
class StaticTimers;
class Corotator;
template <class Scalar> class GenFullSquareMatrix;
typedef GenFullSquareMatrix<double> FullSquareMatrix;
template <class Scalar> class GenFSFullMatrix;
typedef GenFSFullMatrix<double> FSFullMatrix;

extern SfemNonInpc<class Scalar, class VecType> * sfem_noninpc;
extern SfemInpc<class Scalar, class VecType> * sfem_inpc;



template<class T, class VectorType, class SolverType>
class SingleDomainPostProcessor {
 protected:
    Domain *domain;
    T *bcx;
    StaticTimers *times;
    SolverType *solver;
    GenSparseMatrix<T> *kuc, *kcc;
  public:
    SingleDomainPostProcessor(Domain *_domain, T *_bcx, 
                             StaticTimers *_times = NULL, SolverType *_solver = NULL,
                              GenSparseMatrix<T> *_kuc = NULL, GenSparseMatrix<T> *_kcc = NULL)
     : domain(_domain), bcx(_bcx), times(_times), solver(_solver), kuc(_kuc), kcc(_kcc) {}

    void staticOutput(VectorType &solution, VectorType &rhs, bool printTimers = true, int ndflag =0);
    void staticOutput(GeomState& gs, double time = 0.0);
    void staticOutput(GeomState *gs, double time = 0.0);
#ifdef STRUCTOPT
    void staticOutput(VectorType &solution, VectorType &rhs, double time);
#endif
    void getStressStrain(VectorType &sol, int fileNumber, int stressIndex, double time, int printFlag) 
                     { domain->getStressStrain(sol,bcx,fileNumber,stressIndex,time,printFlag); }
    void setsizeSfemStress(int fileNumber) { domain->setsizeSfemStress(fileNumber); }
    int getsizeSfemStress() { return domain->getsizeSfemStress(); }
    T* getSfemStress(int fileNumber) { T* ds = 0; return domain->getSfemStress(fileNumber, ds); }
    void updateSfemStress(T* str, int fileNumber) { domain->updateSfemStress(str, fileNumber); }

    void setSolver(SolverType *s) { solver = s; }
};


template <class T, class VectorType, class SolverType>
class SingleDomainStatic 
{
 protected:
    Domain *domain;
    T *bcx;
    DComplex *bcxC;
    GenSparseMatrix<T> *kuc, *kcc;
    AllOps<T> allOps;
    AllSensitivities<T> allSens;
    SolverType *solver;
    StaticTimers *times; 

    // members for pre-stress linear static problem
    FullSquareMatrix *kelArray;
    Corotator **allCorot;
    GeomState *geomState;

    FSFullMatrix *X;     // pre-calculated projector
    double *Rmem;        // global rigid body modes (numdof X 6)
    int numR;            // number of rigid body modes

 public:
    SingleDomainStatic<T,VectorType,SolverType>(Domain *d) { domain = d; Rmem = 0; numR = 0; }
    int solVecInfo();
    int solVecInfo(int i);
    virtual void getRHS(VectorType &);
    void getRHSinpc(VectorType &);
    void setIWaveDir(int _i) { domain->iWaveDir = _i; }
    void getFreqSweepRHS(VectorType *rhs, VectorType **sol_prev, int iRHS);
    virtual void getRHS(VectorType &,double,double);
    void pade(VectorType *sol, VectorType **sol_prev, double *h, double x) { };
    void preProcessSA();
    void postProcessSA(GenVector<T> &sol);
    virtual void preProcess();
    void rebuildSolver()
      { clean(); preProcess(); }
    void scaleDisp(VectorType &);
    void scaleInvDisp(VectorType &);
    void scaleDisp(VectorType &, double alpha);
    void forceContinuity(VectorType &) {}
    void forceAssemble(VectorType &) {}
    void clean();
    SolverType *getSolver();
    AllOps<T> *getAllOps() { return &allOps; }
    T *getbc() { return bcx; }
    void project(VectorType &f);
    void projector_prep(Rbm *rbms);
    void eigmode_projector_prep();
    SingleDomainPostProcessor<T,VectorType,SolverType> *getPostProcessor();
    StaticTimers *getStaticTimers() { return times; }
    void assignRandMat() {domain->assignRandMat(); }
    void retrieveElemset() {domain->retrieveElemset();}
    AllSensitivities<T> *getAllSensitivities() { return &allSens; }
#ifdef STRUCTOPT
    void preoptProcess();
    void reBuild(); 
    void getPseudoLoad(VectorType &, VectorType &);
    FullSquareMatrix *getkelArray() { return kelArray;} 
#endif
};

#ifdef _TEMPLATE_FIX_
#include <Problems.d/StaticDescr.C>
#endif

#endif
