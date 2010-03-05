#ifndef _NLMATPROBDESC_
#define _NLMATPROBDESC_

#include <Element.d/NLElement.h>
#include <Driver.d/Domain.h>
#include <Driver.d/DynamProbType.h>
#include <Problems.d/StaticDescr.h>   
#include <Driver.d/NLState.h>

class NLMatProbDesc {
   Domain &domain;
   CoordSet &ndset;
   DofSetArray *dsa;
   ConstrainedDSA *c_dsa;
   Solver *solver;
   SparseMatrix *kuc;
   int numNLEle;
   MatNLElement **nlElem;
   int *intStateOffset;
   int numIntDofs;
   int numElemStates;
   double initRes;
   int outIndex;
   int iterNum;
// ------------------------------------- 
    double *bcx;        // displacement prescribed values
    double *vcx;        // velocity prescribed values
    FILE   *res;        // file to store residuals
    int     totIter;     // counter of iterations
    int    *dofTypeArray; //
 
    SparseMatrix *M;    // Mass matrix
    SparseMatrix *Muc;
    SparseMatrix *Mcc;
    Vector *vld;
    Vector *temp;
    Vector *pf;
    Vector *boundAcc;
    Vector *boundVel;
    Vector *boundDsp;
    Vector *boundInertiaF;
    Vector *unconstrInertiaF;
    
    Vector *totalRes;
 
    FullSquareMatrix *kelArray; // array of element stiffness matrices
    FullSquareMatrix *melArray; // array of element mass matrices
 
    PrevFrc *prevFrc;   // previous Aeroelastic force at time step t(n-1)
    double t0;          // initial time
    double totalTime;   // total time
    double dt;          // time step size
    double delta;       // half time step size
    double remainder;   // remaining time after max steps
    int maxStep;        // maximum number of time steps
 
    double tolerance;   // convergence criteria tolerance
    double firstRes;    // first iteration residual norm
    double firstDv;     // first iteration displacement increment norm
    double firstEng;    // first iteration energy increment norm
    double firstForceNorm;
    double firstMomenNorm;
    double externalForceNorm; //current external force norm
 
    int numSystems;     // number of linear systems solved
 
    StaticTimers *times; // Timing information class
    ControlInterface *userSupFunc;
    ControlLawInfo *claw;
 
    void extract(Vector& d_n, Vector& v_n, Vector& a_n,
                 double* ctrdis, double* ctrvel, double* ctracc);
 
    void addCtrl(Vector& externalForce, double *controlForce);
    void addUserForce(Vector & externalForce, double *userForce);   
public:
   NLMatProbDesc(Domain *d);
   void processLastOutput();
   void preProcess();
   int reBuild(int,int,NLState &);
   void staticOutput(NLState *, double, Vector &f, Vector &glRes);
   void printTimers() {}
   void updatePrescribedDisplacement(NLState *, double = 1.0);
   void updatePrescribedDisplacement(NLState *, double, double);
   void getRHS(Vector &, NLState * = 0);
   double getStiffAndForce(NLState &, Vector &, Vector &intrnForce,
		   Vector &glRes);
   Solver *getSolver();
   int elemVecInfo() { return numIntDofs; }
   int solVecInfo() { return numIntDofs; }
   int sysVecInfo() { return dsa->size(); }
   int getMaxit();
   double getMaxLambda();
   int checkConvergence(int, double, double);
   double getDeltaLambda0();
   NLState *createGeomState();
   NLMatProbDesc *getPostProcessor() { return this; }
   
   void updateStates(Vector *internStates, Vector *disp, Vector *prescDisp,
		   Vector *du, Vector *prescDu = 0);

   // Functions for dynamics:
   void computeTimeInfo(); // maybe we need to do something??
   void getConstForce(Vector &gravityForce);
   int  getInitState(Vector& d, Vector& v, Vector& a, Vector &v_p);
   void readRestartFile (Vector &, Vector &, Vector &, Vector &, NLState &);
   double getDelta()  { return delta; }
   double getDt()     { return dt;    }
   void getInitialTime (int &, double &);
   void getExternalForce (Vector &, Vector &, int, double &, NLState *, Vector &, Vector &);
   void dynamCommToFluid (NLState *, NLState *, Vector &, Vector &, Vector &, Vector &,int, int, int) {}
   void dynamOutput (NLState *, Vector &, Vector &, double, int, Vector &, Vector &, Vector &);
   int  getMaxStep() { return maxStep; }
   double getStiffAndForce (NLState &, Vector &, Vector &, double = -1);
   void reBuild (NLState &, int = 0);
   void formRHSpredictor (Vector &, Vector &, Vector &, NLState &, double);
   int checkConvergence (int, double, Vector &, Vector &, double );
   void update (Vector &);
   double integrate(NLState &sn, NLState &snp, Vector &du, Vector &res,
                      Vector &elemForce, Vector &totalRes);
   void printNode (int);
   void get_inc_displacement (Vector &, NLState &);
   double formRHScorrector (Vector &, Vector &, Vector &, Vector &);
   void addInertialTerm(Vector &, Vector &);
   //void midpoint_step_update (Vector &, double &, NLState &);
   void printTimers (double &);
   NLState *copyGeomState(NLState *);
   void setBC(double *userDefineDisplacement, double *userDefineVel); 
   //void updatePrescribedDisplacement(NLState *geomState);  
//   SDDynamPostProcessor *getPostProcessor();      
   double getTolerance()  { return tolerance*firstRes; }
   int getNumStages() { return 1; }
   void setIteration(int i) { iterNum = i; }
   void initNewton() { /* do nothing */ }
   void updateMpcRhs(NLState&) { /* not implemented */ }
   void updateContactConditions(NLState*) { /* not implemented */ }
   void zeroMpcForces() { /* not implemented */ }
   void addMpcForces(Vector &) { /* not implemented */ }
   double norm(Vector &v) { return sqrt(v*v); }
   void deleteContactConditions() { /* not implemented */ }
   void updateSurfaces(NLState*,int) { /* not implemented */ }

   int getAeroAlg() { return domain.solInfo().aeroFlag; }
   int getThermoeFlag() { return domain.solInfo().thermoeFlag; }
   int getThermohFlag() { return domain.solInfo().thermohFlag; }
   int getAeroheatFlag() { return domain.solInfo().aeroheatFlag; }
   void formRHSinitializer(Vector &, Vector &, Vector &, NLState &, Vector &) { cerr << "NLMatProbDesc::formRHSinitializer is not implemented\n"; }
   bool linesearch();
   double getEnergy(double, Vector&, NLState*) { cerr << "NLMatProbDesc::getEnergy is not implemented\n"; }
   void getNewmarkParameters(double &beta, double &gamma, double &alphaf, double &alpham) { cerr << "NLModalDescr::getNewmarkParameters is not implemented\n"; }

};

inline
NLState::NLState(const NLState &state) :
  probDesc(state.probDesc), internalStates(state.internalStates), 
  disp(state.disp), prescDisp(state.prescDisp)
{
}

inline
NLState &NLState::operator=(const NLState &state) {
  internalStates = state.internalStates;
  disp = state.disp;
  prescDisp = state.prescDisp;
  return *this;
}

#endif
