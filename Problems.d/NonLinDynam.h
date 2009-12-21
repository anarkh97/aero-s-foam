#ifndef _NON_LIN_DYNAM_H_
#define _NON_LIN_DYNAM_H_

#include <stdio.h>
#include <stdlib.h>

#include <Driver.d/Domain.h>
#include <Driver.d/DynamProbType.h>

class Domain;
class Rbm;
template <class Scalar> class GenSparseMatrix;
typedef GenSparseMatrix<double> SparseMatrix;
template <class Scalar> class GenSolver;
typedef GenSolver<double> Solver;
template <class Scalar> class GenFullSquareMatrix;
typedef GenFullSquareMatrix<double> FullSquareMatrix;
class StaticTimers;
class GeomState;
class Corotator;
class ControlInterface;
class SDDynamPostProcessor;
template <class Scalar> class GenVector;
typedef GenVector<double> Vector;

#include <Problems.d/StaticDescr.h>

class NLDynamPostProcessor
{
public:
  virtual ~NLDynamPostProcessor() {}
  virtual void dynamCommToFluid(GeomState* geomState, GeomState* bkGeomState, 
                                Vector& velocity, Vector& bkVelocity,
                                Vector& vp, Vector& bkVp, int step, int parity,
                                int aeroAlg) = 0;
  virtual void dynamOutput(GeomState *, GenVector<double> &, GenVector<double> &, double, int, GenVector<double> &, GenVector<double> &,
                           GenVector<double> &) const = 0;
};

// Virtual methods to allow derived class PitaNonLinDynamic in Pita.d/PitaNonLinDynam.d
class NonLinDynamic : public NLDynamPostProcessor {

  protected:
    Domain *domain;
    double *bcx;	// displacement prescribed values
    double *vcx;        // velocity prescribed values
    Solver *solver;
    SparseMatrix *spm;
    Solver *prec;
    SparseMatrix *spp;
    FILE   *res;        // file to store residuals
    int     totIter;     // counter of iterations
    /*int    *dofTypeArray; */ 

    int *clawDofs;      // array containing cdsa dof numbers for usdd nodes

    SparseMatrix *M;    // Mass matrix
    SparseMatrix *C;    // Damping matrix
    SparseMatrix *kuc;
    Corotator **allCorot;
    Vector localTemp;

    FullSquareMatrix *kelArray; // array of element stiffness matrices
    FullSquareMatrix *celArray; // array of element damping matrices
    FullSquareMatrix *melArray; // array of element mass matrices

    PrevFrc *prevFrc;   // previous Aeroelastic force at time step t(n-1)
    double t0;          // initial time
    double totalTime;   // total time
    double dt;          // time step size
    double delta;       // half time step size
    int maxStep;        // maximum number of time steps

    double tolerance;   // convergence criteria tolerance
    double firstRes;    // first iteration residual norm
    double secondRes;    // second iteration residual norm
    double firstDv;     // first iteration displacement increment norm
    double firstEng;    // first iteration energy increment norm
    double firstForceNorm;
    double firstMomenNorm;
    double externalForceNorm; //current external force norm

    int numSystems;     // number of linear systems solved

    StaticTimers *times; // Timing information class

    ControlInterface *userSupFunc;
    ControlLawInfo *claw;

    void extractControlData(Vector& d_n, Vector& v_n, Vector& a_n,
                            double* ctrdis, double* ctrvel, double* ctracc);
    void extractControlDisp(GeomState *, double *);

    void addCtrl(Vector& externalForce, double *controlForce);
    void addUserForce(Vector & externalForce, double *userForce);

    FSFullMatrix *X;    // pre-calculated projector
    double *Rmem;        // global rigid body modes (numdof X 6)
    int numR;            // number of rigid body modes
  
 private:
    virtual void buildOps(AllOps<double> &, double, double, double, Rbm*);

 public:
    // Constructor
    NonLinDynamic(Domain *d);
    virtual ~NonLinDynamic();

    SDDynamPostProcessor * getPostProcessor();
    virtual const NLDynamPostProcessor & defaultPostProcessor() const;
    void getInitialTime(int &initTimeIndex, double &initTime);
    void readRestartFile(Vector &d_n, Vector &v_n, Vector &a_n,
                         Vector &v_p, GeomState &geomState);
    void setBC(double *userDefineDisplacement, double *userDefineVel);

    int  getInitState(Vector& d, Vector& v, Vector& a, Vector &v_p);
    void updateUserSuppliedFunction(Vector& d_n, Vector& v_n, Vector &a_n, Vector &v_p, double initialTime);
    void updatePrescribedDisplacement(GeomState *geomState);

    int  solVecInfo() { return domain->numUncon(); } // number of unconstrained dof
    int  sysVecInfo() { return domain->numdof(); }
    int  elemVecInfo();

    double getTolerance(){ return (tolerance*firstRes); }

    void trProject(Vector &f);
    void projector_prep(Rbm *R, SparseMatrix *M);

    void   computeTimeInfo();

    double getDelta() const     { return delta;     }
    double getDt() const        { return dt;        }
    int    getMaxStep() const   { return maxStep;   }
    double getTotalTime() const { return totalTime; }
    int    getMaxit();
    double getDeltaLambda();

    void getConstForce(Vector & gravityForce);
    void getExternalForce(Vector & externalForce, 
                          Vector & gravityForce, int tIndex, double time,
                          GeomState* geomState, Vector& elementInternalForce, 
                          Vector &aeroF);

    double formRHScorrector(Vector& inc_displac, Vector &velocity, Vector& acceleration,
                            Vector &residual, Vector &rhs);
    double formRHScorrector(Vector& inc_displac, Vector &velocity, Vector& acceleration,
                            Vector &residual, Vector &rhs, double localDelta);

    void formRHSpredictor(Vector &velocity, Vector &acceleration, Vector &residual, Vector &rhs, GeomState &, double mid = 0.0);
    void formRHSpredictor(Vector &velocity, Vector &acceleration, Vector &residual, Vector &rhs, GeomState &, double mid, double localDelta);

    void formRHSinitializer(Vector &fext, Vector &velocity, Vector &elementInternalForce, GeomState &geomState, Vector &rhs);

    void preProcess();

    virtual void openResidualFile();

    void processLastOutput();
    Solver *getSolver();

    const SparseMatrix* getMassMatrix() { return M; }

    GeomState* createGeomState();
    GeomState* copyGeomState(GeomState* geomState);

    int getNumStages();
    int checkConvergence(int iter, double rN, Vector& residual, Vector& dv, double time);

    // getStiffAndForce forms element stiffness matrices and
    // returns the residual force = external - internal forces
    double getStiffAndForce(GeomState& geomState, Vector& residual, 
                          Vector& elementInternalForce, double midtime=-1);


    // reBuild assembles new dynamic stiffness matrix
    void reBuild(GeomState& geomState, int iter = 0);
    void reBuild(GeomState& geomState, int iter, double localDelta);

    void printTimers(double timeLoop);
    virtual void dynamCommToFluid(GeomState* geomState, GeomState* bkGeomState,
                          Vector& velocity, Vector& bkVelocity,
                          Vector& vp, Vector& bkVp, int step, int parity,
                          int aeroAlg);
    virtual void dynamOutput(GeomState* geomState, Vector& velocity, Vector &vp,
                     double time, int timestep, Vector& force, Vector &aeroF, Vector &acceleration) const;
    void initNewton() { /* do nothing */ }

    int  aeroPreProcess(Vector& d_n, Vector& v_n, Vector& a_n, Vector& vp);
    void a5TimeLoopCheck(int& parity, double& t, double dt);
    void a5StatusRevise(int parity, SysState<Vector>& curState, SysState<Vector>& bkState);
    int getAeroAlg() { return domain->solInfo().aeroFlag; }
    int getThermoeFlag() { return domain->solInfo().thermoeFlag; }
    void getNewmarkParameters(double &beta, double &gamma,
                              double &alphaf, double &alpham);

};

inline double
NonLinDynamic::formRHScorrector(Vector &inc_displacement, Vector &velocity, Vector& acceleration,
                                Vector &residual,         Vector &rhs)
{
  return formRHScorrector(inc_displacement, velocity, acceleration, residual, rhs, delta);
}

inline void
NonLinDynamic::formRHSpredictor(Vector &velocity, Vector &acceleration, Vector &residual, Vector &rhs, GeomState &geomState, double mid)
{
  formRHSpredictor(velocity, acceleration, residual, rhs, geomState, mid, delta);
}

inline void
NonLinDynamic::reBuild(GeomState& geomState, int iter)
{
  reBuild(geomState, iter, delta);
}

inline const NLDynamPostProcessor &
NonLinDynamic::defaultPostProcessor() const
{
  return *this;
}

class DummyNLDynamPostProcessor {
public:
  virtual void dynamCommToFluid(GeomState* geomState, GeomState* bkGeomState,
                          Vector& velocity, Vector& bkVelocity,
                          Vector& vp, Vector& bkVp, int step, int parity,
                          int aeroAlg){}
  virtual void dynamOutput(GeomState *, GenVector<double> &, GenVector<double> &, double, int, GenVector<double> &, GenVector<double> &,
                           GenVector<double> &) const {}
};

#endif
