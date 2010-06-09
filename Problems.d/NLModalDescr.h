#ifndef _NL_MODAL_DESCR_H_
#define _NL_MODAL_DESCR_H_

#include <Problems.d/ModalBase.h>
#include <Problems.d/ModalGeomState.h>
#include <Math.d/DBSparseMatrix.h>
#include <Math.d/matrix.h>

class NLModalOpSolver{
/* wrapper class for a FullMatrix
   passed to NLDynamSolver as template parameter OpSolver
*/
private:

  FullM mat;

public:

  NLModalOpSolver(){}
  void reSolve(Vector &rhs);

  friend class NLModalDescr;

};

//------------------------------------------------------------------------------
//******************************************************************************
//------------------------------------------------------------------------------

class NLModalDescr : public ModalBase{
/* problem descriptor for non-linear dynamics problem
   solves eq of motion in a body-fixed reference frame
*/

private:

  Vector fullFrot;

  Vector **modesPr;  // rotational modes of the deformed structure
  double J[3][3];    // mass moment of inertia tensor

  int numFilters;
  int *rbmIdx;         // indicies of the rbms to be filtered
  FullM *invPsiTMPsi;  // inverse of Psi^T.M.Psi used for the rbm filter
  FullM *MPsi;         // M.Psi; Psi is the matrix of rbms to be filtered
  

  double ****Bp;     // the various B quantities are coupling coefficients
  double  ***Bh;     // note to self: Bp is B'; Bh, Bhat; Bt, Btilde
  double  ***Bt;

  double ***Bq;      // temporary variables to store intermediate products
  double ***Bo;
  double Bqq[3][3];
  double Bhq[3][3];
  double **Bho;
  double **Btq;
  double **Bqo;
  double **Bto;
  Vector q_modesPr[3];

// TODO: remove these data members ---------------
  double **dxhat;       // dxhat and dflexhat are derivatives wrt theta
  double ***dflexhat;
  double **flexhat;     // for each flex mode, global dofs involved in a constr

  double **dxhatnp1;    // same as above but evaluated at n+1 time step
  double ***dflexhatnp1;
  double **flexhatnp1;
//------------------------------------------------

  double *constr;       // the constraints evaluated at n+1 state

// new stuff 040902 ------------------------------
  double **dhdp;   // derivative of the constraints wrt generlaized translation
  double **dhdth;  //  "  "  "  wrt generalized rotation
  double **dhdq;   //  "  "  "  wrt modal deformation
  double ***d2h;   // 2nd derivative wrt modal deformation and rotation

  double **dhdp_np1;   // same derivatives as above, but evaluated
  double **dhdth_np1;  //   at the n+1  time step
  double **dhdq_np1;
// -----------------------------------------------


  double c22a[3][3];    // temporary variables to store intermediate sums
  double c22b[3];       // c22b is used in the rhs and in NLModalOpSolver.mat
  double mcoef, ccoef;  // coefficients for the contribution of M, C and K in
  double kcoef, hcoef;  //   NLModalOpSolver.mat; hcoef is scaling factor for constraints

  NLModalOpSolver *solver;
  double tFinal;
  double dt;
  double delta;
  int    maxStep;
  double remainingTime;
  double tolerance;     // tolerance for convergence criteria
  double firstRes;      // norm of the residual of the first iteration of
                        //   a given time step
public:

  NLModalDescr(){}
  NLModalDescr(Domain *d);

  void rotateVector(Vector& globalF, Vector& rotF, double glR[3][3]);
  void projectForceNL(Vector& fullF, Vector& modalF, ModalGeomState* mgs);
  void filterForce(Vector& fullF);
  void expandFlex(Vector& modalV, Vector& fullV);
  void calcTmps(ModalGeomState &mgs);

  void preProcess();
  void processLastOutput();

  void buildRBMFilter(DBSparseMatrix* massMat);

  void computeTimeInfo();
  NLModalOpSolver* getSolver(){ return solver; }

  int solVecInfo(){ return (numRBM + numFlex + numConstr); }
  int sysVecInfo(){ return domain->numdof(); }
  int elemVecInfo(){ return domain->maxNumDOF(); }

  void getConstForce(Vector &constF);
  void getExternalForce(Vector &extF, Vector &gravF, int tIndex,
    double time, ModalGeomState* mgs, Vector &elemIntF, Vector& aeroF);

  int getInitState(Vector &dsp, Vector &vel, Vector &acc, Vector &vel_p);

  ModalGeomState* createGeomState(){ return new ModalGeomState(numRBM, numFlex, numConstr); }
  ModalGeomState* copyGeomState(ModalGeomState* mgs)
    { return new ModalGeomState(*mgs); }

  void readRestartFile(Vector &dsp, Vector &vel, Vector &acc, Vector &vel_p,
    ModalGeomState &mgs);
  void updatePrescribedDisplacement(ModalGeomState *mgs);
  
  int    getMaxit()  { return domain->solInfo().getNLInfo().maxiter; }
  int    getMaxStep(){ return maxStep; }
  double getDt()     { return dt; }
  double getDelta()  { return delta; }

  void getInitialTime(int &initTimeIndex, double &initTime){
    initTimeIndex = domain->solInfo().initialTimeIndex;
    initTime      = domain->solInfo().initialTime;
  }

  void initIterState(ModalGeomState &mgs);
  double getStiffAndForce(ModalGeomState &mgs,
    Vector &res, double midtime = -1);
  void reBuild(ModalGeomState &mgs, int iter = 0){ solver->mat.factor(); }

  void evalRHS(Vector &res, Vector &rhs, ModalGeomState &mgs);
  void formRHSinitializer(Vector &, Vector &, Vector &, ModalGeomState &, Vector &) { cerr << "NLModalDescr::formRHSinitializer is not implemented\n"; }
  void formRHSpredictor(Vector &res, Vector &rhs, ModalGeomState &mgs);
  double formRHScorrector(Vector& res, Vector& rhs, ModalGeomState& mgs);

  int checkConvergence(int iter, double normRes, Vector &residual,
    Vector &dvec, double time);

  double getTolerance(){ return (tolerance*firstRes); }

  void dynamCommToFluid(ModalGeomState* mgs, ModalGeomState* bkMgs, 
                        Vector& vel, Vector& bkVel,
                        Vector& vel_p, Vector& bkVel_p,
                        int tIndex, int parity,
                        int aeroAlg){}
  void dynamOutput(ModalGeomState* mgs, Vector& vel, Vector& vel_p,
    double time, int tIndex, Vector& extF, Vector &aeroF, Vector &acc);
  void printTimers(double timeLoop){ /* leave blank for now */ }

  int getNumStages(){ return 1; }

  void dRdTheta(double R[3][3], double dR[3][3][3]);
  void evalConstraintsAndDerivatives(ModalGeomState &mgs);
  
  void test2(ModalGeomState* = 0);
  void test(ModalGeomState* = 0);
  void printCoefs();
  double getResidualNorm(Vector &rhs) { return rhs.norm(); }

  int getAeroAlg() { return domain->solInfo().aeroFlag; }
  int getThermoeFlag() { return domain->solInfo().thermoeFlag; }
  int getThermohFlag() { return domain->solInfo().thermohFlag; }
  int getAeroheatFlag() { return domain->solInfo().aeroheatFlag; }

  void getNewmarkParameters(double &beta, double &gamma, double &alphaf, double &alpham) { cerr << "NLModalDescr::getNewmarkParameters is not implemented\n"; }
};

#endif
