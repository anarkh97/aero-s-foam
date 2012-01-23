#ifndef _MD_NL_DYNAMIC_H_
#define _MD_NL_DYNAMIC_H_

#include <map>
#include <vector>
#include <cstddef>

class Domain;
template <class Scalar> class GenDecDomain;
typedef GenDecDomain<double> DecDomain;
template <class Scalar> class GenParallelSolver;
typedef GenParallelSolver<double> ParallelSolver;
template <class Scalar> class GenSparseMatrix;
typedef GenSparseMatrix<double> SparseMatrix;
template <class Scalar> class GenFullSquareMatrix;
typedef GenFullSquareMatrix<double> FullSquareMatrix;
class Corotator;
class StaticTimers;
class ControlInterface;
class ControlLawInfo;
template <class Scalar> class GenDistrVector;
typedef GenDistrVector<double> DistrVector;
class DistrGeomState;
template <class Scalar> class GenMultiDomainPostProcessor;
typedef GenMultiDomainPostProcessor<double> MultiDomainPostProcessor;
template <class Scalar> class GenSubDOp;
typedef GenSubDOp<double> SubDOp;
template <class Scalar> class GenMDDynamMat;
typedef GenMDDynamMat<double> MDDynamMat;
class DistrInfo;
class DistFlExchanger;

// Multiple Domain Nonlinear Dynamic problem descriptor

class MDNLDynamic 
{
    Domain *domain;
    DecDomain  *decDomain;
    ParallelSolver *solver;
    int     totIter;             // counter of iterations

    MDDynamMat *allOps;
    SubDOp *M;
    SubDOp *C;
    SubDOp *Kuc;
    Corotator ***allCorot;       // element corotators per subdomain
    DistrVector *localTemp;

    FullSquareMatrix **kelArray; // element stiffness matrices per subdomain
    FullSquareMatrix **melArray; // element mass matrices per subdomain
    FullSquareMatrix **celArray; // element damping matrices per subdomain

    double t0;                   // initial time
    double totalTime;		 // total simulation time
    double dt;			 // time step size
    double delta;		 // half time step size
    double remainder;		 // remaining time after max steps
    int    maxStep;		 // maximum number of time steps

    double tolerance;		 // newton convergence criteria tolerance
    double firstRes;		 // first iteration residual norm
    double secondRes;            // second iteration residual norm
    double firstDv;		 // first iteration displacement increment norm
    double firstEng;             // first iteration energy increment norm
    double externalForceNorm;    // current external force norm

    int numSystems;              // number of linear systems solved

    ControlInterface *userSupFunc;
    ControlLawInfo *claw;
    int **clawDofs;

    StaticTimers *times;	 // timing information

    // user defined displacements and velocities
    double **usrDefDisps;
    double **usrDefVels;

    // aero data
    DistFlExchanger *distFlExchanger;
    DistrVector *prevFrc;
    int prevIndex;
    double prevTime;
    DistrVector *aeroForce;
    DistrVector *nodalTemps;

    std::map<std::pair<int,int>, double> *mu; // lagrange multipliers for the contact surfaces
    std::vector<double> *lambda; // lagrange multipliers for all other constraints

 public:

    // Constructor
    MDNLDynamic(Domain *d);
    virtual ~MDNLDynamic();

    MultiDomainPostProcessor *getPostProcessor();
    void getInitialTime(int &initTimeIndex, double &initTime);
    void readRestartFile(DistrVector &d_n, DistrVector &v_n, DistrVector &a_n,
                         DistrVector &v_p, DistrGeomState &geomState); 

    int getInitState(DistrVector &d, DistrVector& v, DistrVector &a, DistrVector &v_p);
    void updatePrescribedDisplacement(DistrGeomState *geomState);

    DistrInfo& solVecInfo();
    DistrInfo& elemVecInfo();
    DistrInfo& sysVecInfo();

    double getTolerance() { return (tolerance*firstRes); }

    void computeTimeInfo();

    double getDelta() const     { return delta; }
    double getDt() const        { return dt;    }
    double getLastDt() const    { return remainder; }
    int    getMaxStep() const   { return maxStep;   }
    double getTotalTime() const { return totalTime; }
    int    getMaxit();
    double getDeltaLambda();

    void getConstForce(DistrVector &gravityForce);
    void getExternalForce(DistrVector &externalForce, DistrVector &constantForce,
                          int tIndex, double time, DistrGeomState *geomState, 
                          DistrVector &elementInternalForce, DistrVector &aeroF);

    double formRHScorrector(DistrVector& inc_displacement, DistrVector& velocity, DistrVector& acceleration,
                           DistrVector& residual, DistrVector& rhs);
    double formRHScorrector(DistrVector& inc_displacement, DistrVector& velocity, DistrVector& acceleration,
                           DistrVector& residual, DistrVector& rhs, double localDelta);
    
    void formRHSpredictor(DistrVector& velocity, DistrVector& acceleration, DistrVector& residual, 
                          DistrVector& rhs, DistrGeomState &, double mid = 0.0);
    void formRHSpredictor(DistrVector& velocity, DistrVector& acceleration, DistrVector& residual,
                          DistrVector& rhs, DistrGeomState &, double mid, double localDelta);

    void formRHSinitializer(DistrVector &fext, DistrVector &velocity, DistrVector &elementInternalForce,
                            DistrGeomState &geomState, DistrVector &rhs, DistrGeomState *refState = NULL);

    virtual void preProcess();

    void processLastOutput();
    ParallelSolver *getSolver();

    DistrGeomState* createGeomState();
    DistrGeomState* copyGeomState(DistrGeomState* geomState);

    int getNumStages();
    int checkConvergence(int iter, double rN, DistrVector& residual, DistrVector& dv, double time);

    void updateStates(DistrGeomState *refState, DistrGeomState& geomState);

    // getStiffAndForce forms element stiffness matrices and
    // returns the residual force = external - internal forces
    double getStiffAndForce(DistrGeomState& geomState, DistrVector& residual,
                            DistrVector& elementInternalForce, double midtime=-1, DistrGeomState *refState = NULL);

    // reBuild assembles new dynamic stiffness matrix
    void reBuild(DistrGeomState& geomState, int iter = 0);
    void reBuild(DistrGeomState& geomState, int iter, double localDelta);

    void printTimers(double timeLoop);
    void dynamOutput(DistrGeomState* geomState, DistrVector& velocity, DistrVector &vp,
                     double time, int timestep, DistrVector& force, DistrVector &aeroF, DistrVector &acceleration);
    void dynamCommToFluid(DistrGeomState* geomState, DistrGeomState* bkGeomState,
                          DistrVector& velocity, DistrVector& bkVelocity,
                          DistrVector& vp, DistrVector& bkVp, int step, int parity,
                          int aeroAlg);

    double getResidualNorm(DistrVector &vec);

    int getAeroAlg();
    int getThermoeFlag();
    int getThermohFlag();
    int getAeroheatFlag();
    void getNewmarkParameters(double &beta, double &gamma,
                              double &alphaf, double &alpham);
  protected:
    Domain *getDomain() { return domain; }
    DecDomain *getDecDomain() { return decDomain; }

  private:
    void makeSubDofs(int isub);
    void makeSubCorotators(int isub);
    void makeSubElementArrays(int isub);
    void subGetExternalForce(int isub, DistrVector& f, DistrVector& constantForce, double time);
    void subGetStiffAndForce(int isub, DistrGeomState &geomState,
                             DistrVector &res, DistrVector &elemIntForce, double t, DistrGeomState *refState);
    void subUpdatePrescribedDisplacement(int isub, DistrGeomState& geomState);
    void addConstraintForces(int isub, DistrVector& rhs);
    void getConstraintMultipliers(int isub);
    void subUpdateGeomStateUSDD(int isub, DistrGeomState &geomState, double *userDefineDisp);
    void makeSubClawDofs(int isub);
    void subKucTransposeMultSubtractClaw(int iSub, DistrVector& residual, double *userDefineDisp);
    void subExtractControlDisp(int isub, DistrGeomState &geomState, double *ctrdsp);
    int aeroPreProcess(DistrVector &, DistrVector &, DistrVector &, DistrVector &);
    void thermoePreProcess();
    void thermohPreProcess(DistrVector &);
    void aeroheatPreProcess(DistrVector &, DistrVector &, DistrVector &);
    void subDynamCommToFluid(int isub, DistrVector& v, DistrGeomState* distrGeomState,
                             DistrGeomState* bkDistrGeomState, int parity, int aeroAlg);
    void subDynamCommToFluidAeroheat(int isub, DistrVector& v, DistrGeomState* distrGeomState);
    void updateConstraintTerms(DistrGeomState* geomState, double t);
    void subUpdateStates(int isub, DistrGeomState *refState, DistrGeomState *geomState);
    void subReadRestartFile(int i, DistrVector &d_n, DistrVector &v_n, DistrVector &a_n,
                            DistrVector &v_p, DistrGeomState &geomState);
    void subWriteRestartFile(int i, double &t, int &index, DistrVector &vel_n, DistrGeomState &geomState);

    virtual bool factorWhenBuilding() const;
};

inline double
MDNLDynamic::formRHScorrector(DistrVector& inc_displacement, DistrVector& velocity, DistrVector& acceleration,
                              DistrVector& residual, DistrVector& rhs)
{
  return formRHScorrector(inc_displacement, velocity, acceleration, residual, rhs, delta);
}

inline void 
MDNLDynamic::formRHSpredictor(DistrVector& velocity, DistrVector& acceleration, DistrVector& residual,
                              DistrVector& rhs, DistrGeomState &geomState, double mid)
{
  formRHSpredictor(velocity, acceleration, residual, rhs, geomState, mid, delta);
}

inline void
MDNLDynamic::reBuild(DistrGeomState& geomState, int iter)
{
  reBuild(geomState, iter, delta);
}

#endif
