#ifndef _MD_DYNAM_DESCR_H_
#define _MD_DYNAM_DESCR_H_

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
   Rbm* rigidBodyModes;

   GenMDDynamMat() { dynMat = 0; Msolver = 0; K = 0; C = 0; Cuc = 0; M = 0; Muc = 0; Mcc = 0; Kuc = 0; C_deriv = 0; Cuc_deriv = 0; rigidBodyModes = 0; };
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
    DistrVector *nodalTemps;

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
    void setNodalTemps(DistrVector*);
    void dynamOutput(int, MDDynamMat &, DistrVector &, DistrVector *aeroF, SysState<DistrVector> &);
    double getKineticEnergy(DistrVector & vel, SubDOp * gMass) { return 0.0; }
};

class MultiDomainDynam 
{
protected:
    DecDomain *decDomain;
    Domain *domain;

private:
    CuCSparse **cucs;
    StaticTimers *times;

    // control law data
    ControlInterface *userSupFunc;
    ControlLawInfo *claw;

    GenFullSquareMatrix<double> **kelArray;
    Corotator ***allCorot;
    DistrGeomState *geomState;
    DistrVector *dprev;

    MDDynamMat *dynMat;
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

    const DistrInfo &solVecInfo();
    DistrInfo &bcInfo();

    int getTimeIntegration();
    int getFilterFlag();
    int* boundary();
    double* boundaryValue();
    Domain* getDomain();
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

    void trProject(DistrVector &);
    void project(DistrVector &);
    void getRayleighCoef(double& alpha);

    void addPrescContrib(SubDOp*, SubDOp*, DistrVector&, DistrVector&, 
                         DistrVector&, DistrVector&, double t);

    SubDOp* getpK(MDDynamMat* dynOps);
    SubDOp* getpM(MDDynamMat* dynOps);
    SubDOp* getpC(MDDynamMat* dynOps);

    // Central Difference only related subroutines
    void computeStabilityTimeStep(double, MDDynamMat&);

    // Mode Decomposition parameters and subroutines
    int getModeDecompFlag();
    void modeDecompPreProcess(SparseMatrix* M);
    void modeDecomp(double t, int tIndex, DistrVector& d_n);

    void getInternalForce(DistrVector &d, DistrVector &f, double t);

    // Aeroelastic problems related subroutines
    void computeTimeInfo();
    int aeroPreProcess(DistrVector &, DistrVector &, DistrVector &, 
		       DistrVector &); 
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
    void subGetInternalForce(int isub, DistrVector &res, double t);
    void subGetKtimesU(int isub, DistrVector &d, DistrVector &f);
    void makeSubCorotators(int isub);
    void makeSubElementArrays(int isub);
    void initSubPrescribedDisplacement(int isub);
    void subUpdateGeomStateUSDD(int isub, double *userDefineDisp);
    void subExplicitUpdate(int isub, DistrVector &d);
};

#endif
