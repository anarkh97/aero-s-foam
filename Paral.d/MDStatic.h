#ifndef _MD_STATIC_DESCR_H_
#define _MD_STATIC_DESCR_H_
#include <Driver.d/DecDomain.h>

template <class Scalar> class GenFetiSolver;
class StaticTimers;
template <class Scalar> class GenDistrVector;
class DistrGeomState;
class DistrInfo;

template <class Scalar> 
class GenMultiDomainPostProcessor 
{
 protected:
    GenDecDomain<Scalar> *decDomain;
    GenFetiSolver<Scalar> *solver;
    StaticTimers *times;

  public:
    GenMultiDomainPostProcessor(GenDecDomain<Scalar> *d, GenFetiSolver<Scalar> *s, 
                                StaticTimers* _times=0)
       { decDomain = d; solver = s; times = _times; }
    void staticOutput(GenDistrVector<Scalar> &, GenDistrVector<Scalar> &, bool printTimers = true, int ndflag=0);
    void staticOutput(DistrGeomState *u, double lambda = 1.0);
//    void staticOutput_convertdispstress(GenDistrVector<Scalar> &,  GenDistrVector<Scalar> &, int fileNumber,
//                   int stressIndex, double time = 0) { cerr << "GenMultiDomainPostProcessor::staticOutput_convertdispstress not implemented" << endl; }
    void getStressStrain(GenDistrVector<Scalar> &, int fileNumber,
                   int stressIndex, double time, int pflag) { cerr << "GenMultiDomainPostProcessor::getStressStrain not implemented" << endl; }
    void setsizeSfemStress(int fileNumber) { cerr << "GenMultiDomainPostProcessor::setsizeSfemStress(int fileNumber) not implemented" << endl; }
    int getsizeSfemStress() { cerr << "GenMultiDomainPostProcessor::getsizeSfemStress() not implemented" << endl; return 0; }
    Scalar* getSfemStress(int fileNumber) {cerr << "GenMultiDomainPostProcessor::getSfemStress() not implemented" << endl; return 0;}
    void updateSfemStress(Scalar* str, int fileNumber) { cerr << "GenMultiDomainPostProcessor::updateSfemStress() not implemented" << endl;}
    void getMemoryK(int iSub, long *memory);
    void getMemoryPrec(int iSub, long *memory);
};

template<class Scalar>
class GenMultiDomainStatic 
{
 protected:
    Domain *domain;
    GenDecDomain<Scalar> *decDomain;
    GenFetiSolver<Scalar> *solver;
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
    void getRHSinpc(GenDistrVector<Scalar>  &) {cerr << "GenMultiDomainStatic::getRHSinpc not implemented" << endl ; }
    void preProcess();
    void assignRandMat() {cerr << "GenMultiDomainStatic::assignRandMat() not implemented" << endl ; }
    void retrieveElemset() {cerr << "GenMultiDomainStatic::retrieveElemset() not implemented" << endl ; }
    void rebuildSolver();
    void scaleDisp(GenDistrVector<Scalar> &);
    void scaleInvDisp(GenDistrVector<Scalar> &);

    void clean();
    void setIWaveDir(int _i); // FETI-H
    void getFreqSweepRHS(GenDistrVector<Scalar> *rhs, GenDistrVector<Scalar> **sol_prev, int iRHS);
    void getRHS(GenDistrVector<Scalar> &,double,double);
    void pade(GenDistrVector<Scalar> *sol, GenDistrVector<Scalar> **sol_prev, double *h, double x);
    GenMDDynamMat<Scalar> *getAllOps() { return &allOps; }
    GenFetiSolver<Scalar> *getSolver();
    GenMultiDomainPostProcessor<Scalar> *getPostProcessor();
    StaticTimers *getStaticTimers() { return times; }
    void project(GenDistrVector<Scalar> &); // { fprintf(stderr,"Paral.d/MDStatic.h: project not implemented here"); }
 private:
    void eigmode_projector_prep();

};

typedef GenMultiDomainStatic<double> MultiDomainStatic;
typedef GenMultiDomainPostProcessor<double> MultiDomainPostProcessor;

#ifdef _TEMPLATE_FIX_
#include <Paral.d/MDStatic.C>
#endif

#endif
