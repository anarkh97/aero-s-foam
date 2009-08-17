#ifndef _MD_DYNAM_DESCR_H_
#define _MD_DYNAM_DESCR_H_

#include <Driver.d/DynamProbType.h>
#include <Solvers.d/ParallelSolver.h>

template <class Scalar> class GenVector;
typedef GenVector<double> Vector;
template <class Scalar> class GenSparseMatrix;
typedef GenSparseMatrix<double> SparseMatrix;
template <class Scalar> class GenCuCSparse;
typedef GenCuCSparse<double> CuCSparse;
class StaticTimers;
class DistFlExchanger;
class ControlInterface;

template <class Scalar> class GenDecDomain;
typedef GenDecDomain<double> DecDomain;
template <class Scalar> class GenFetiSolver;
typedef GenFetiSolver<double> FetiSolver;
template <class Scalar> class GenDistrVector;
typedef GenDistrVector<double> DistrVector;
template <class Scalar> class GenSubDOp;
typedef GenSubDOp<double> SubDOp;
class Domain;

//template <class VecType> class SysState;

template<class Scalar>
class GenMDDynamMat {
 public:
   union {
     GenParallelSolver<Scalar> *dynMat;
     GenParallelSolver<Scalar> *sysSolver;
   };
   GenParallelSolver<Scalar> *Msolver;
   GenSubDOp<Scalar> *K;
   GenSubDOp<Scalar> *C;
   GenSubDOp<Scalar> *Cuc;
   GenSubDOp<Scalar> *M;
   GenSubDOp<Scalar> *Muc;
   GenSubDOp<Scalar> *Mcc;
   GenSubDOp<Scalar> *Kuc;
   GenSubDOp<Scalar> **C_deriv;
   GenSubDOp<Scalar> **Cuc_deriv;

   GenMDDynamMat() { dynMat = 0; Msolver = 0; K = 0; C = 0; Cuc = 0; M = 0; Muc = 0; Mcc = 0; Kuc = 0; C_deriv = 0; Cuc_deriv = 0; };
/*
   ~MDDynamMat() { 
     if(dynMat) delete dynMat; 
     if(K) delete K; 
     if(C) delete C; 
     if(Cuc) delete Cuc; 
     if(M) delete M;
     if(Muc) delete Muc;
     if(Mcc) delete Mcc;
     if(Kuc) Kuc->partialClean();
   }
*/
};

typedef GenMDDynamMat<double> MDDynamMat;

class MultiDomDynPostProcessor 
{
    DistFlExchanger *distFlExchanger;
 
    // user defined displacements and velocities
    double **usrDefDisps;
    double **usrDefVels;
    GenDecDomain<double> *decDomain;
    StaticTimers *times;

  public:
    MultiDomDynPostProcessor(DecDomain *d, StaticTimers* _times = 0) {
      decDomain = d;
      times = _times;
    }
    MultiDomDynPostProcessor(DecDomain *d, 
		DistFlExchanger *_distFlExchanger, StaticTimers* _times = 0) {
      decDomain = d;
      distFlExchanger = _distFlExchanger;
      times = _times;
    }
    void setPostProcessor(DistFlExchanger *);
    void setUserDefs(double **, double **);
    void dynamOutput(int, MDDynamMat &, DistrVector &, DistrVector *aeroF, SysState<DistrVector> &);
    double getKineticEnergy(DistrVector & vel, SubDOp * gMass) { return 0.0; }
};

template <class Scalar>
class MultiDomainDynam 
{
    DecDomain *decDomain;
    Domain *domain;
    CuCSparse **cucs;
    FetiSolver *solver;
    //MDDynamMat dynMat;
    StaticTimers *times;

    // control law data
    ControlInterface *userSupFunc;
    ControlLawInfo *claw;

    FullSquareMatrix **kelArray;
    Corotator ***allCorot;
    DistrGeomState *geomState;
    DistrVector *dprev;

  public:
    MultiDomainDynam(Domain *d);
    ~MultiDomainDynam();
    MDDynamMat * buildOps(double, double, double);
    MultiDomDynPostProcessor *getPostProcessor();

    const DistrInfo &solVecInfo();
    DistrInfo &bcInfo();

    int getTimeIntegration();
    int getFilterFlag();
    int *boundary() { return 0; }

    double * boundaryValue() { return 0; }
	
    Domain *getDomain() { return domain; }
    void getTimes(double &dt, double &t);
    void getNewMarkParameters(double &beta, double &gamma,
                              double &alphaf, double &alpham);
    void getQuasiStaticParameters(double &maxVel, double &delta);
    void getInitState(SysState<DistrVector> &);
    void getInitialTime(int &tIndex, double &initialTime);
    double getInitialForceNorm();
    void getSteadyStateParam(int &steadyFlag, int &steadyMin, int &steadMax,
                             double &steadyTol); 
    void getConstForce(DistrVector &);
    void getContactForce(DistrVector &, DistrVector &ctc_f);
    void computeExtForce2(SysState<DistrVector> &, DistrVector &, 
                          DistrVector &, int tIndex, double t,
                          DistrVector * aero_f=0,
                          double gamma=0.5, double alphaf=0.5);

    void getRHS(DistrVector &);
    void preProcess();
    void processLastOutput();
    void printTimers(MDDynamMat *, double);

    void trProject(DistrVector &) { fprintf(stderr,"Paral.d/MDDynam.h: trProject not implemented here"); }
    void project(DistrVector &) { fprintf(stderr,"Paral.d/MDDynam.h: project not implemented here"); }
    void computeStabilityTimeStep(double, MDDynamMat &);
    void getRayleighCoef(double &alpha) { alpha = domain->solInfo().alphaDamp; }

    void addPrescContrib(SubDOp*, SubDOp*, DistrVector&, DistrVector&, 
                         DistrVector&, DistrVector&, double t );

    SubDOp * getpK(MDDynamMat * dynOps) { return dynOps->K; }
    SubDOp * getpM(MDDynamMat * dynOps) { return dynOps->M; }
    SubDOp * getpC(MDDynamMat * dynOps) { return dynOps->C; }

    void getInternalForce(DistrVector &d, DistrVector &f, double t);
    // these can be private
    void subGetInternalForce(int isub, DistrVector &res);
    void subGetKtimesU(int isub, DistrVector &d, DistrVector &f);
    void makeSubCorotators(int isub);
    void makeSubElementArrays(int isub);
    void initSubPrescribedDisplacement(int isub);
    void subUpdateGeomStateUSDD(int isub, double *userDefineDisp);

    // Aeroelastic problems related subroutines
    // these are not defined yet. AEROELASTIC is not defined to be used
    // with FETI solver yet.
  private:
    DistFlExchanger *distFlExchanger; 
    MultiDomDynPostProcessor *mddPostPro;
    
    // user defined displacements and velocities
    double **usrDefDisps;
    double **usrDefVels;

    // Previous Force
    DistrVector *prevFrc;
    int prevIndex;
    double prevTime;

    // Aero Force
    DistrVector *aeroForce;
 
  public:
    void computeTimeInfo() { cerr << "MultiDomainDynam::computeTimeInfo() is not implemented\n"; };
    int aeroPreProcess(DistrVector &, DistrVector &, DistrVector &, 
		       DistrVector &); 
    void a5TimeLoopCheck(int &, double &, double) { cerr << "MultiDomainDynam::a5TimeLoopCheck is not implemented\n"; };
    void a5StatusRevise(int, SysState<DistrVector> &, SysState<DistrVector> &) { cerr << "MultiDomainDynam::a5StatusRevise is not implemented\n"; };
    int getAeroAlg() { return domain->solInfo().aeroFlag; }
    int getThermoeFlag() { return domain->solInfo().thermoeFlag; }
    void aeroSend(double time, DistrVector& d, DistrVector& v, DistrVector& a, DistrVector& v_p) 
        { cerr << "MultiDomainDynam::aeroSend() is not implemented\n"; };

    int cmdCom(int cmdFlag); 
   
    // data member access functions
    double **getUserDefDisps()  { return usrDefDisps; }
    double **getUserDefVels()   { return usrDefVels; }

    void thermoePreProcess(DistrVector &, DistrVector &, DistrVector &) {}
    void modeDecompPreProcess(SparseMatrix *M) {}
    void modeDecomp(double t, int tIndex, DistrVector& d_n) {}
    int getModeDecompFlag() { return 0; }
};
#ifdef _TEMPLATE_FIX_
  #include <Paral.d/MDDynamTem.C>
//  #include <Paral.d/MDDynam.C>
#endif

#endif
