#ifndef _MD_DYNAM_DESCR_H_
#define _MD_DYNAM_DESCR_H_

#include <cstdlib>
#ifdef DISTRIBUTED
#include <Utils.d/DistHelper.h>
#endif

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
template <class Scalar> class GenDistrVector;
typedef GenDistrVector<double> DistrVector;
template <typename Scalar> class GenFullSquareMatrix;
template <class Scalar> class GenSubDOp;
typedef GenSubDOp<double> SubDOp;
class Domain;
class Rbm;

template <typename T> class SysState;
template <typename T> class GenParallelSolver;
template <typename Scalar> struct AllSensitivities;

class ControlLawInfo;
class Corotator;
class DistrInfo;
class DistrGeomState;

template<class Scalar>
class GenMDDynamMat {
 public:
   union {
     GenParallelSolver<Scalar> *dynMat;
     GenParallelSolver<Scalar> *sysSolver;
   };
   GenSubDOp<Scalar> *spMat;

   GenParallelSolver<Scalar> *prec;       // preconditioner
   GenSubDOp<Scalar> *spp;

   GenParallelSolver<Scalar> *Msolver;
   GenSubDOp<Scalar> *K;
   GenSubDOp<Scalar> *Kuc;
   GenSubDOp<Scalar> *Kcc;
   GenSubDOp<Scalar> *C;
   GenSubDOp<Scalar> *Cuc;
   GenSubDOp<Scalar> *Ccc;
   GenSubDOp<Scalar> *M;
   GenSubDOp<Scalar> *Muc;
   GenSubDOp<Scalar> *Mcc;
   GenSubDOp<Scalar> **C_deriv;
   GenSubDOp<Scalar> **Cuc_deriv;
   Rbm* rigidBodyModes;

   GenMDDynamMat() { dynMat = 0; spMat = 0; prec = 0; spp = 0; Msolver = 0; K = 0; Kuc = 0; Kcc = 0; C = 0; Cuc = 0; Ccc = 0; M = 0; Muc = 0; Mcc = 0; 
                     C_deriv = 0; Cuc_deriv = 0; rigidBodyModes = 0; };
};

typedef GenMDDynamMat<double> MDDynamMat;

class MultiDomDynPostProcessor 
{
  protected:
    DistFlExchanger *distFlExchanger;
 
    // user defined displacements and velocities
    double **usrDefDisps;
    double **usrDefVels;
    GenDecDomain<double> *decDomain;
    StaticTimers *times;
    DistrVector *nodalTemps;
    DistrGeomState *geomState;
    Corotator ***allCorot;

  public:
    MultiDomDynPostProcessor(DecDomain *d, StaticTimers* _times, DistrGeomState *_geomState = 0,
                             Corotator ***_allCorot = 0) {
      decDomain = d;
      times = _times;
      geomState = _geomState;
      allCorot = _allCorot;
    }
    MultiDomDynPostProcessor(DecDomain *d, DistFlExchanger *_distFlExchanger, StaticTimers* _times,
                             DistrGeomState *_geomState = 0, Corotator ***_allCorot = 0) {
      decDomain = d;
      distFlExchanger = _distFlExchanger;
      times = _times;
      geomState = _geomState;
      allCorot = _allCorot;
    }
    void setPostProcessor(DistFlExchanger *);
    void setUserDefs(double **, double **);
    void setNodalTemps(DistrVector*);
    void dynamOutput(int, double, MDDynamMat &, DistrVector &, DistrVector *aeroF, SysState<DistrVector> &);
    double getKineticEnergy(DistrVector & vel, SubDOp * gMass) { return 0.0; }
};

class MultiDomainDynam 
{
protected:
    DecDomain *decDomain;
    Domain *domain;
    AllSensitivities<double> *allSens;

private:
    CuCSparse **cucs;
    StaticTimers *times;

    // control law data
    ControlInterface *userSupFunc;
    ControlLawInfo *claw;

protected:
    GenFullSquareMatrix<double> **kelArray;
    GenFullSquareMatrix<double> **melArray;
    Corotator ***allCorot;
    DistrGeomState *geomState;
    MDDynamMat *dynMat;

private:
    MultiDomDynPostProcessor *mddPostPro;

    // user defined displacements and velocities
    double **usrDefDisps;
    double **usrDefVels;

    // aero data
    DistFlExchanger *distFlExchanger;
    DistrVector *aeroForce;
    DistrVector *prevFrc;
    int prevIndex;
    double prevTime;
    DistrVector *prevFrcBackup;
    int prevIndexBackup;
    double prevTimeBackup;

    // thermoe/thermoh data
    DistrVector* nodalTemps;

  public:
    MultiDomainDynam(Domain *d);
    virtual ~MultiDomainDynam();
    MDDynamMat * buildOps(double, double, double);
    MultiDomDynPostProcessor *getPostProcessor();

    const DistrInfo &solVecInfo() const;
    DistrInfo &bcInfo();

    int getTimeIntegration();
    int getFilterFlag();
    int* boundary();
    double* boundaryValue();
    Domain* getDomain();
    AllSensitivities<double> *getAllSensitivities() { return allSens; }
    void getTimes(double &dt, double &t);
    void getNewMarkParameters(double &beta, double &gamma,
                              double &alphaf, double &alpham);
    void getQuasiStaticParameters(double &maxVel, double &delta);
    void getInitState(SysState<DistrVector> &);
    void printFullNorm(DistrVector &){};
    void getInitialTime(int &tIndex, double &initialTime);
    double getInitialForceNorm();
    void getSteadyStateParam(int &steadyFlag, int &steadyMin, int &steadMax,
                             double &steadyTol); 
    void getConstForce(DistrVector &);
    void addConstForceSensitivity(DistrVector &);
    void getContactForce(DistrVector &d_n, DistrVector &dinc, DistrVector &ctc_f, double t_n_p, double dt, double dt_old);
    void computeExtForce2(SysState<DistrVector> &, DistrVector &, 
                          DistrVector &, int tIndex, double t,
                          DistrVector * aero_f=0,
                          double gamma=0.5, double alphaf=0.5);
    void getAeroelasticForceSensitivity(int t_index, double t, DistrVector * aero_f=0, double gamma=0.5, double alphaf=0.5) {
      filePrint(stderr," ... MultiDomainDynam::getAeroelasticForceSensitivity\n");  exit(-1);
    }

    void getRHS(DistrVector &);
    void preProcess();
    void preProcessSA() {  filePrint(stderr," ... MultiDomainDynam::preProcessSA is not implemented\n");  exit(-1);  }
    void postProcessSA(DistrVector &sol) {  filePrint(stderr," ... MultiDomainDynam::postProcessSA is not implemented\n");  exit(-1);  }
    void processLastOutput();
    void printTimers(MDDynamMat *, double);

    void trProject(DistrVector &);
    void project(DistrVector &);
    void getRayleighCoef(double& alpha);

    void addPrescContrib(SubDOp*, SubDOp*, DistrVector&, DistrVector&, 
                         DistrVector&, DistrVector&, double tm, double tf);

    SubDOp* getpK(MDDynamMat* dynOps);
    SubDOp* getpM(MDDynamMat* dynOps);
    SubDOp* getpC(MDDynamMat* dynOps);

    // Central Difference only related subroutines
    virtual void computeStabilityTimeStep(double&, MDDynamMat&);

    void updateState(double dt_n_h, DistrVector& v_n_h, DistrVector& d_n);

    // Mode Decomposition parameters and subroutines
    int getModeDecompFlag();
    void modeDecompPreProcess(SparseMatrix* M);
    void modeDecomp(double t, int tIndex, DistrVector& d_n);

    void getInternalForce(DistrVector &d, DistrVector &f, double t, int tIndex);

    // Aeroelastic problems related subroutines
    void computeTimeInfo();
    int aeroPreProcess(DistrVector &, DistrVector &, DistrVector &, DistrVector &); 
    int aeroSensitivityPreProcess(DistrVector &, DistrVector &, DistrVector &, DistrVector &); 
    int sendDisplacements(DistrVector &, DistrVector &, DistrVector &, DistrVector &); 
    int cmdCom(int cmdFlag);
    int getAeroAlg();
    void aeroSend(double time, DistrVector& d, DistrVector& v, DistrVector& a, DistrVector& v_p);
    void a5TimeLoopCheck(int&, double&, double);
    void a5StatusRevise(int, SysState<DistrVector>&, SysState<DistrVector>&);

    // Thermoelastic problems related subroutines
    void thermoePreProcess(DistrVector&, DistrVector&, DistrVector&);
    int getThermoeFlag();
    void thermohPreProcess(DistrVector&, DistrVector&, DistrVector&);
    int getThermohFlag();

    // Aeroheat
    void aeroHeatPreProcess(DistrVector&, DistrVector&, DistrVector&);
    int getAeroheatFlag();
   
  private:
    void subGetInternalForce(int isub, DistrVector &res, double &t, int &tIndex);
    void subGetKtimesU(int isub, DistrVector &d, DistrVector &f);
    void makeSubCorotators(int isub);
    void makeSubElementArrays(int isub);
    void initSubPrescribedDisplacement(int isub);
    void subUpdateGeomStateUSDD(int isub, double *userDefineDisp, DistrGeomState *geomState,
                                double *userDefineVel, double *userDefineAcc);
    void subExplicitUpdate(int isub, DistrVector &d, DistrGeomState *geomState);
};

#endif
