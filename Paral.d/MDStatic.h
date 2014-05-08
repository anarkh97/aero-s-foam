#ifndef _MD_STATIC_DESCR_H_
#define _MD_STATIC_DESCR_H_
#include <Driver.d/DecDomain.h>

template <class Scalar> class GenParallelSolver;
class StaticTimers;
template <class Scalar> class GenDistrVector;
class DistrGeomState;
class DistrInfo;
template <class Scalar> struct AllSensitivities;

template <class Scalar> 
class GenMultiDomainPostProcessor 
{
 protected:
    GenDecDomain<Scalar> *decDomain;
    GenParallelSolver<Scalar> *solver;
    StaticTimers *times;

  public:
    GenMultiDomainPostProcessor(GenDecDomain<Scalar> *d, GenParallelSolver<Scalar> *s, 
                                StaticTimers* _times=0)
       { decDomain = d; solver = s; times = _times; }
    void staticOutput(GenDistrVector<Scalar> &, GenDistrVector<Scalar> &, bool printTimers = true, int ndflag=0);
    void staticOutput(DistrGeomState *u, double lambda = 1.0);
    void getStressStrain(GenDistrVector<Scalar> &, int fileNumber,
                   int stressIndex, double time, int pflag);
    void setsizeSfemStress(int fileNumber);
    int getsizeSfemStress();
    Scalar* getSfemStress(int fileNumber);
    void updateSfemStress(Scalar* str, int fileNumber);
    void getMemoryK(int iSub, long *memory);
    void getMemoryPrec(int iSub, long *memory);
    void setSolver(GenParallelSolver<Scalar> *s) { solver = s; }
};

template<class Scalar>
class GenMultiDomainStatic 
{
 protected:
    Domain *domain;
    GenDecDomain<Scalar> *decDomain;
    GenParallelSolver<Scalar> *solver;
    StaticTimers *times;
    GenMDDynamMat<Scalar> allOps;
//    SubDOp *K;
//    SubDOp *M;
//    SubDOp **C_deriv;
    int numR;
    GenDistrVector<Scalar> **Rmem;
    GenFSFullMatrix<Scalar> *RtRinverse;
 public:
    GenMultiDomainStatic() : decDomain(0), solver(0), times(0) {}
    explicit GenMultiDomainStatic(Domain *d);
    ~GenMultiDomainStatic() { delete decDomain; delete times; }  // solver deleted in StaticSolver

    DistrInfo &solVecInfo();
    DistrInfo &solVecInfo(int i);
    void getRHS(GenDistrVector<Scalar> &);
    void getRHSinpc(GenDistrVector<Scalar> &);
    void preProcessSA();
    void postProcessSA(GenDistrVector<Scalar> &);
    void preProcess();
    void assignRandMat();
    void retrieveElemset();
    void rebuildSolver();
    void scaleDisp(GenDistrVector<Scalar> &);
    void scaleInvDisp(GenDistrVector<Scalar> &);
    void scaleDisp(GenDistrVector<Scalar> &, double alpha);
    void forceContinuity(GenDistrVector<Scalar> &);
    void forceAssemble(GenDistrVector<Scalar> &);

    void clean();
    void setIWaveDir(int _i); // FETI-H
    void getFreqSweepRHS(GenDistrVector<Scalar> *rhs, GenDistrVector<Scalar> **sol_prev, int iRHS);
    void getRHS(GenDistrVector<Scalar> &,double,double);
    void pade(GenDistrVector<Scalar> *sol, GenDistrVector<Scalar> **sol_prev, double *h, double x);
    GenMDDynamMat<Scalar> *getAllOps() { return &allOps; }
    GenParallelSolver<Scalar> *getSolver();
    GenMultiDomainPostProcessor<Scalar> *getPostProcessor();
    StaticTimers *getStaticTimers() { return times; }
    void project(GenDistrVector<Scalar> &);
    AllSensitivities<Scalar> *getAllSensitivities() { cerr << "GenMultiDomainStatic::getAllSensitivities() not implemented" << endl; return 0; }
 private:
    void eigmode_projector_prep();
    void subGetRHS(int isub, GenDistrVector<Scalar>& rhs);
    void makeSubdomainStaticLoadGalPr(int iSub, GenDistrVector<Scalar> &f, GenDistrVector<Scalar> &tmp, double *o);
    void subPade(int iSub, GenDistrVector<Scalar> *sol, GenDistrVector<Scalar> **u, double *h, double x);
    void multM(int iSub, GenDistrVector<Scalar> *rhs, GenDistrVector<Scalar> **u, int k);
    void multMCoupled1(int iSub, GenDistrVector<Scalar> *rhs, GenDistrVector<Scalar> **u, int k);
    void multMCoupled2(int iSub, GenDistrVector<Scalar> *rhs);
};

typedef GenMultiDomainStatic<double> MultiDomainStatic;
typedef GenMultiDomainPostProcessor<double> MultiDomainPostProcessor;

#ifdef _TEMPLATE_FIX_
#include <Paral.d/MDStatic.C>
#endif

#endif
