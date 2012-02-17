#include <Driver.d/NLMat.h>

#include <Math.d/FullSquareMatrix.h>
#include <Math.d/Skyline.d/SkyMatrix.h>
#include <Math.d/SparseMatrix.h>
#include <Math.d/DBSparseMatrix.h>
#include <Math.d/NBSparseMatrix.h>
#include <Math.d/CuCSparse.h>

#include <Math.d/Vector.h>
#include <Control.d/ControlInterface.h>
#include <Corotational.d/GeomState.h>
#include <Corotational.d/Corotator.h>
#include <Hetero.d/FlExchange.h>

extern int verboseFlag;

NLMatProbDesc::NLMatProbDesc(Domain *d):
 domain(*d), ndset(d->getNodes()) {
  outIndex = 0;
  numSystems = 0;
  temp = 0;
  pf   = 0;
  dofTypeArray = 0;
  tolerance = domain.solInfo().getNLInfo().tolRes;
}

void
NLMatProbDesc::processLastOutput()  {

  OutputInfo *oinfo = geoSource->getOutputInfo();
  for (int iOut = 0; iOut < geoSource->getNumOutInfo(); iOut++)
    oinfo[iOut].interval = 1;
}

void
NLMatProbDesc::preProcess()
{
 // Allocate space for the Static Timers
 times = new StaticTimers;

  domain.preProcessing();
  domain.make_constrainedDSA();
  dsa = domain.getDSA();
  c_dsa = domain.getCDSA();
  domain.makeAllDOFs();
  numIntDofs = c_dsa->size();
  Elemset *allElems = domain.getEset();
  int nEle = allElems->size();
  
  nlElem = new MatNLElement *[nEle];
  int iEle;
  // Grab all the non linear elements
  numNLEle = 0;
  for(iEle = 0; iEle < nEle; ++iEle) {
    if((*allElems)[iEle]) {
       MatNLElement *thisElem = 
            dynamic_cast<MatNLElement *>( (*allElems)[iEle] );
       if(thisElem != 0) {
         nlElem[numNLEle++] = thisElem;
       }
    }
  }
 if(nEle != numNLEle)
    fprintf(stderr, "Not all elements are non linear. This may not be supported\n");
 intStateOffset = new int[numNLEle];

 numElemStates = 0;
 for(iEle = 0; iEle < numNLEle; ++iEle) {
    intStateOffset[iEle] = numElemStates;
    numElemStates += nlElem[iEle]->numStates();
 }

 // Look if there is a user supplied routine for control
 claw = geoSource->getControlLaw();

 userSupFunc = geoSource->getUserSuppliedFunction();
  
 if(domain.probType() == SolverInfo::MatNonLinStatic) {
/* PJSA 
    domain.getSolverAndKuc<double>(solver, kuc);
*/
    AllOps<double> allOps;

    domain.buildOps<double>(allOps, 1.0, 0.0, 0.0, (Rbm *) NULL, (FullSquareMatrix *) NULL,
                            (FullSquareMatrix *) NULL, true);

    solver = allOps.sysSolver;
    spm = allOps.spm;
    
    domain.createKelArray(kelArray);
 }
 else {
   AllOps<double> allOps;

   allOps.M = domain.constructDBSparseMatrix<double>();
   allOps.Kuc = domain.constructCuCSparse<double>();
   allOps.Muc = domain.constructCuCSparse<double>();
   allOps.Mcc = domain.constructCCSparse<double>();

   double Kcoef = 0.0; // PJSA changed to 0.0 from 1.0 
   double Mcoef = 1.0;
   double Ccoef = 0.0;

   domain.buildOps<double>(allOps, Kcoef, Mcoef, Ccoef, (Rbm *) NULL, (FullSquareMatrix *) NULL,
                           (FullSquareMatrix *) NULL, true);

   M      = allOps.M;
   Mcc    = allOps.Mcc;
   Muc    = allOps.Muc;
   solver = allOps.sysSolver;
   spm    = allOps.spm;
   kuc = allOps.Kuc;

   boundAcc = new Vector(Mcc->numRow(), 0.0);
   boundVel = new Vector(Mcc->numRow(), 0.0);
   boundDsp = new Vector(Mcc->numRow(), 0.0);
   boundInertiaF = new Vector(Mcc->numRow(), 0.0);
   unconstrInertiaF = new Vector(M->numRow(), 0.0);
   
   // melArray is only needed for dynamics
   domain.createKelArray(kelArray, melArray);
 }
 totalRes = new Vector(sysVecInfo());
}

Solver *
NLMatProbDesc::getSolver()
{
 return solver;
}

void
NLMatProbDesc::getRHS(Vector &force, NLState *state)
{
 // buildRHS with a zero kuc will not add the non homogeneous forces
 domain.buildRHSForce<double>(force);
}

void
NLMatProbDesc::updateStates(Vector *internStates, Vector *disp,
                   Vector *prescDisp,
		   Vector *du, Vector *prescDu)
{
  int iState = 0;
  double *states = internStates->data();
  int iEle;
  for(iEle = 0; iEle < numNLEle; ++iEle) {
    int dofs[128]; // XML Check this is not too small
    int nd = nlElem[iEle]->numDofs();
    Node nodes[32];
    int ndNum[32];
    double un[128];
    double unp[128];
    if(nd > 128)
      throw "Error in maximum number of dofs";
    int nnd = nlElem[iEle]->numNodes();
    if(nnd > 32)
      throw "Error in maximum number of nodes";
    nlElem[iEle]->nodes(ndNum);
    for(int iNode = 0; iNode < nnd; ++iNode)
      nodes[iNode] = *ndset[ndNum[iNode]];
    nlElem[iEle]->dofs(*dsa, dofs);
    int iDof;
    for(iDof = 0; iDof < nd; ++iDof) {
      int df = c_dsa->getRCN(dofs[iDof]);
      if(df >= 0) {
	un[iDof] = (*disp)[df];
	unp[iDof] = (*disp)[df] + (*du)[df];
      } else {
	int fdf = c_dsa->invRCN(dofs[iDof]);
	un[iDof] = (*prescDisp)[fdf];
	unp[iDof] = (*prescDisp)[fdf] + ((prescDu == 0) ? 0 : (*prescDu)[fdf]);
      }
    }
    nlElem[iEle]->updateStates(nodes, states+iState, un, unp);
    iState += nlElem[iEle]->numStates();
  }

 //disp->print("Un\n");
}

NLState *
NLMatProbDesc::createGeomState()
{
  NLState *res = new NLState(this,
		  numElemStates, numIntDofs, dsa->size()-numIntDofs);
  res->disp.zero();
  res->prescDisp.zero();
  int iState = 0;
  double *states = res->internalStates.data();
  int iEle;
  for(iEle = 0; iEle < numNLEle; ++iEle) {
    nlElem[iEle]->initStates(states+iState);
    iState += nlElem[iEle]->numStates();
  }
 
  return res;
}

#define NEWWAY

double
NLMatProbDesc::integrate(NLState &sn, NLState &snp, Vector &du, Vector &resid,
		         Vector &intrnForce, Vector &glRes)
{
  int iState = 0;
  double *internStateN = sn.internalStates.data();
  double *internStateNp = snp.internalStates.data();
  int iEle, iDof;
  snp.disp = sn.disp + du;

  //fprintf(stderr, "Integrating\n");
  //state.prescDisp.print("Presc Disp\n");
  //state.disp.print("Current disp\n");
  // zero the force
  intrnForce.zero();
  SparseMatrix *skm = dynamic_cast<SparseMatrix *> (solver);
  if(skm == 0) {
    fprintf(stderr,
	"Error in NLMatProbDesc::getStiffAndForce\n");
  }
#ifndef NEWWAY
  skm->zeroAll();
  kuc->zeroAll();
#endif
  glRes.zero();
  for(iDof = 0; iDof < dsa->size(); ++iDof) {
    int df = c_dsa->getRCN(iDof);
    if(df >= 0)
      glRes[iDof] = resid[df];
  }
  for(iEle = 0; iEle < numNLEle; ++iEle) {
    int iDof, iNode;
    int dofs[128]; // XML Check this is not too small
    int nd = nlElem[iEle]->numDofs();
    Node nodes[32];
    int ndNum[32];
    double un[128];
    double unp[128];
    double force[128];
    //double locK[128*128];
    if(nd > 128)
      throw "Error in maximum number of dofs";
    int nnd = nlElem[iEle]->numNodes();
    if(nnd > 32)
      throw "Error in maximum number of nodes";
    nlElem[iEle]->nodes(ndNum);
    for(iNode = 0; iNode < nnd; ++iNode)
      nodes[iNode] = *ndset[ndNum[iNode]];
    nlElem[iEle]->dofs(*dsa, dofs);
    for(iDof = 0; iDof < nd; ++iDof) {
      int df = c_dsa->getRCN(dofs[iDof]);
      if(df >= 0) {
	un[iDof]  = sn.disp[df];
	unp[iDof] = snp.disp[df];
      } else {
	int fdf = c_dsa->invRCN(dofs[iDof]);
	un[iDof]  = sn.prescDisp[fdf];
	unp[iDof] = snp.prescDisp[fdf];
      }
    }
    //FullSquareMatrix k(nd, locK);
    FullSquareMatrix &kel = kelArray[iEle];
    // nlElem[iEle]->getStiffAndForce(nodes, un, internStates+iState,
    //                                fsm, force);
    nlElem[iEle]->integrate(nodes, un,  internStateN+iState, 
                                   unp, internStateNp+iState,
                            kel, force);
 
    for(iDof = 0; iDof < nd; ++iDof) {
      int df = c_dsa->getRCN(dofs[iDof]);
      if(df >= 0)
	resid[df] += force[iDof];
      glRes[dofs[iDof]] +=  force[iDof];
     }
#ifndef NEWWAY
    skm->add(kel, dofs);
    kuc->add(kel, dofs);
#endif
    iState += nlElem[iEle]->numStates();
  }
  //resid.print("Residual\n");
#ifndef NEWWAY
  solver->factor();
  //fprintf(stderr, "Factored and found %d RBMs\n", solver->numRBM());
#endif
  if(solver->numRBM() > 0) {
     fprintf(stderr,
         "Error: The tangent stiffness matrix has %d singularities\n",
              solver->numRBM());
     exit(-10);
  }

  // pressure using surfacetopo
  domain.addPressureForce(resid);

  return resid.norm();
}

double
NLMatProbDesc::getStiffAndForce(NLState &state, Vector &resid,
		Vector &intrnForce, Vector &glRes)
{
  int iState = 0;
  double *internStates = state.internalStates.data();
  int iEle, iDof;
  
  //state.prescDisp.print("Presc Disp\n");
  //state.disp.print("Current disp\n");
  // zero the force
  intrnForce.zero();
  SparseMatrix *skm = dynamic_cast<SparseMatrix *> (solver);
  if(skm == 0) {
    fprintf(stderr,
	"Error in NLMatProbDesc::getStiffAndForce\n");
  }
  skm->zeroAll();
  kuc->zeroAll();
  glRes.zero();
  for(iDof = 0; iDof < dsa->size(); ++iDof) {
    int df = c_dsa->getRCN(iDof);
    if(df >= 0)
      glRes[iDof] = resid[df];
  }
  for(iEle = 0; iEle < numNLEle; ++iEle) {
    int iDof, iNode;
    int dofs[128]; // XML Check this is not too small
    int nd = nlElem[iEle]->numDofs();
    Node nodes[32];
    int ndNum[32];
    double un[128];
    double force[128];
    double locK[128*128];
    if(nd > 128)
      throw "Error in maximum number of dofs";
    int nnd = nlElem[iEle]->numNodes();
    if(nnd > 32)
      throw "Error in maximum number of nodes";
    nlElem[iEle]->nodes(ndNum);
    for(iNode = 0; iNode < nnd; ++iNode)
      nodes[iNode] = *ndset[ndNum[iNode]];

    nlElem[iEle]->dofs(*dsa, dofs);
    for(iDof = 0; iDof < nd; ++iDof) {
      int df = c_dsa->getRCN(dofs[iDof]);
      if(df >= 0) {
	un[iDof] = state.disp[df];
      } else {
	int fdf = c_dsa->invRCN(dofs[iDof]);
	un[iDof] = state.prescDisp[fdf];
      }
    }
    FullSquareMatrix fsm(nd, locK);
    nlElem[iEle]->getStiffAndForce(nodes, un, internStates+iState,
		    fsm, force);
    for(iDof = 0; iDof < nd; ++iDof) {
      int df = c_dsa->getRCN(dofs[iDof]);
      if(df >= 0)
	resid[df] += force[iDof];
      glRes[dofs[iDof]] +=  force[iDof];
     }
    skm->add(fsm, dofs);
    kuc->add(fsm, dofs);
    iState += nlElem[iEle]->numStates();
  }
  //resid.print("Residual\n");
  solver->factor();
  if(solver->numRBM() > 0) {
     fprintf(stderr,
         "Error: The tangent stiffness matrix has %d singularities\n",
              solver->numRBM());
     exit(-10);
  }
  return resid.norm();
}

void
NLMatProbDesc::updatePrescribedDisplacement(NLState *state, double lambda)
{
  domain.buildPrescDisp(state->prescDisp, lambda);
  if(claw != 0) {
   double *userDefinedDisp = (double *) dbg_alloca(sizeof(double)*claw->numUserDisp);
   double *userDefinedVel  = (double *) dbg_alloca(sizeof(double)*claw->numUserDisp);
   userSupFunc->usd_disp( lambda, userDefinedDisp, userDefinedVel );
   int i;
   for(i=0; i<claw->numUserDisp; ++i) {

     int dsaDof = dsa->locate(claw->userDisp[i].nnum,
                         1 << claw->userDisp[i].dofnum);
     if(dsaDof < 0)
       continue;
     int cdof = c_dsa->invRCN(dsaDof);
     if(cdof >= 0)
      state->prescDisp[cdof] = userDefinedDisp[i];
   }
  }
}

void
NLMatProbDesc::updatePrescribedDisplacement(NLState *state, double t,
		double dt)
{ 
  double coef = 4.0/dt;
  domain.buildPrescDisp(state->prescDisp, t, dt);
  if(claw != 0) {
   double *userDefinedDisp = (double *) dbg_alloca(sizeof(double)*claw->numUserDisp);
   double *userDefinedVel  = (double *) dbg_alloca(sizeof(double)*claw->numUserDisp);
   userSupFunc->usd_disp( t, userDefinedDisp, userDefinedVel );
   int i;
   for(i=0; i<claw->numUserDisp; ++i) {

     int dsaDof = dsa->locate(claw->userDisp[i].nnum,
                         1 << claw->userDisp[i].dofnum);
     if(dsaDof < 0)
       continue;
     int cdof = c_dsa->invRCN(dsaDof);
     if(cdof >= 0)
      state->prescDisp[cdof] = userDefinedDisp[i];
   }
  }

  int size = Mcc->numRow();
  int i;
  for(i = 0; i < size; ++i)
    boundVel[i] = coef*(state->prescDisp[i] - (*boundDsp)[i])
		   -(*boundVel)[i];
  for(i = 0; i < size; ++i) {
    double t = state->prescDisp[i];
    state->prescDisp[i] = 2*state->prescDisp[i]-(*boundDsp)[i];
    (boundAcc)[i] = 0.25*coef*coef*(state->prescDisp[i]-(*boundDsp)[i])
		    -coef*(*boundVel)[i]-(*boundAcc)[i];
    (*boundDsp)[i] = t;
  }

}

int
NLMatProbDesc::checkConvergence(int it, double, double resNorm)
{
  //fprintf(stderr, "Init res %e\n", initRes);
  if(it == 0)
    initRes = resNorm;
  double tolerance = domain.solInfo().getNLInfo().tolRes; 
  fprintf(stderr, "Current: %e Init res: %e tolerance: %e\n", 
            resNorm, initRes, tolerance);
  if(resNorm <= tolerance*initRes) return 1;
  return 0;
}

void
NLMatProbDesc::staticOutput(NLState *state, double t, Vector &f, Vector &totRes, NLState *)
{
  domain.postProcessing<double>(state->disp,state->prescDisp.data(), f, 0, outIndex, t);
  domain.resProcessing(totRes, outIndex, t);
  outIndex++;
}

int
NLMatProbDesc::getMaxit()
{
  return domain.solInfo().getNLInfo().maxiter;
}

double
NLMatProbDesc::getDeltaLambda0()
{
 return domain.solInfo().getNLInfo().dlambda;
}
 
double
NLMatProbDesc::getMaxLambda()
{
 return domain.solInfo().getNLInfo().maxLambda;
}

bool
NLMatProbDesc::linesearch()
{
 return domain.solInfo().getNLInfo().linesearch;
}

NLState::NLState(NLMatProbDesc*pdesc,int nInt,int nfree,int npresc) :
		probDesc(pdesc),
                internalStates(nInt), disp(nfree), prescDisp(npresc)
{
  disp.zero();
  prescDisp.zero();
}

void
NLState::update(NLState &refState, Vector&du)
{
 internalStates = refState.internalStates;
 disp = refState.disp;
 // Vector finalPrescDisp(prescDisp);
 Vector dPrU(prescDisp);
 dPrU -= refState.prescDisp;
 // prescDisp = refState.prescDisp;
 // This should work because updateStates does not touch the presc disp vector
 //du.print("du");
 //dPrU.print("dPrU");
 //internalStates.print("Before");
 probDesc->updateStates(&internalStates,&disp,&refState.prescDisp, &du, &dPrU); 
 disp += du;
 //internalStates.print("After");
 // prescDisp = finalPrescDisp;
}


// FUNCTIONS FOR DYNAMICS
void NLMatProbDesc::computeTimeInfo()
{
  // Time integration information
 
  // Get total time and time step size and store them
  totalTime = domain.solInfo().tmax;
  dt        = domain.solInfo().getTimeStep();
  delta     = 0.5*dt;
 
  // Compute maximum number of steps
  maxStep = (int) ( totalTime/dt );
 
  // Compute time remainder
  remainder = totalTime - maxStep*dt;
 
  // set half time step size in user defined functions
  if(userSupFunc)
    userSupFunc->setDt(delta);
}

void
NLMatProbDesc::getConstForce(Vector& constantForce)
{
 domain.computeConstantForce(constantForce);
/*
 constantForce.zero();
 
 if( domain.gravityFlag()  ) domain.buildGravityForce<double>(constantForce);
 
 // if( domain.pressureFlag() ) domain.buildPressureForce(constantForce);
*/
}

int
NLMatProbDesc::getInitState(Vector& d_n, Vector& v_n, Vector &a_n, Vector &v_p)
{
  // Initialize displacement, velocity and acceleration vectors
 domain.initDispVeloc(d_n, v_n, a_n, v_p);
 
 // For getting a user supplied control function
 if(userSupFunc) {
   double *ctrdisp = (double *) dbg_alloca(sizeof(double)*claw->numSensor);
   double *ctrvel  = (double *) dbg_alloca(sizeof(double)*claw->numSensor);
  
 double *ctracc  = (double *) dbg_alloca(sizeof(double)*claw->numSensor);
   extract( d_n, v_n, a_n, ctrdisp, ctrvel, ctracc );
   userSupFunc->init( ctrdisp, ctrvel, ctracc );
 }
 return -1;
}

void
NLMatProbDesc::readRestartFile(Vector &d_n, Vector &v_n, Vector &a_n,
                               Vector &v_p, NLState &geomState)
{
  d_n.zero();
  v_n.zero();
  a_n.zero();
  v_p.zero();
  // MLXX This need to be updated for material non-linear
  // domain.readRestartFile(d_n, v_n, a_n, v_p, bcx, vcx, geomState);       
}

NLState*
NLMatProbDesc::copyGeomState(NLState* geomState)
{
 return new NLState(*geomState);
}

void
NLMatProbDesc::getInitialTime(int &initTimeIndex, double &initTime)
{
 initTimeIndex = domain.solInfo().initialTimeIndex;
 initTime      = domain.solInfo().initialTime;  
}

void
NLMatProbDesc::getExternalForce(Vector& rhs, Vector& gf, int tIndex, double& t,
                            NLState* geomState, Vector& elemNonConForce, Vector &aeroF)
{
 // ... BUILD THE RHS FORCE (external + gravity + nonhomogeneous)
 times->formRhs -= getTime();
 
 if(domain.pressureFlag()) {
   // allocate Vector pf only one time
   if(pf == 0) pf = new Vector(solVecInfo());
   Vector &pressureForce = *pf;
   pressureForce.zero();
   // MLXX This need to be updated for material non-linear
   // domain.buildPressureForce( pressureForce, geomState);
   gf += pressureForce;
 }
 domain.computeExtForce4(rhs, gf, t, kuc);

/* XXXXX XML Work on this for user-supplied functions
 // Compute any user defined forces
 if(claw > 0) {
   double *userDefineForce=(double *) dbg_alloca(sizeof(double)*claw->numUserForce);
   userSupFunc->usd_forc(t,userDefineForce);
   addUserForce(rhs, userDefineForce);
 }

 // Compute any user defined control forces
 if(userSupFunc) {

   double *ctrdisp = (double *) dbg_alloca(sizeof(double)*claw->numSensor);
   double *ctrvel  = (double *) dbg_alloca(sizeof(double)*claw->numSensor);
   double *ctracc  = (double *) dbg_alloca(sizeof(double)*claw->numSensor);
   double *ctrfrc  = (double *) dbg_alloca(sizeof(double)*claw->numActuator);

   // KHP: we need the state of the control sensors to pass to
   //      the user supplied control function
   // extract(d_n, v_n, a_n, ctrdisp,ctrvel,ctracc);

   userSupFunc->ctrl(ctrdisp,ctrvel,ctracc,ctrfrc);
   addCtrl(rhs, ctrfrc);
 }
*/

}
void
NLMatProbDesc::dynamOutput(NLState* state, Vector& velocity,
                           Vector& vp, double time, int step, 
                           Vector& force, Vector &aeroF, Vector &acceleration)
{  
  //if(solver->numRBM() != 0)
  //  fprintf(stderr,"\n --- Number of RBM %4d             ...\n",
  //         solver->numRBM());
 
  times->output -= getTime();
  
 
  if(step == -1) 
    vld = new Vector(domain.maxNumDOF());

  domain.postProcessing<double>(state->disp,state->prescDisp.data(), force, 0, step+1, time);
 // domain.resProcessing(totRes, outIndex, t);
 times->output += getTime();

}

/*int
NLMatProbDesc::getMaxStep()
{

}
*/

double
NLMatProbDesc::getStiffAndForce(NLState& state, Vector& residual,
                                Vector& elementInternalForce, double midtime)
{
 times->buildStiffAndForce -= getTime();
 
 double *userDefinedDisp = 0;
 double *userDefinedVel  = 0;
 
 if( midtime != -1.0 && claw != 0 ) {
   userDefinedDisp = (double *) dbg_alloca(sizeof(double)*claw->numUserDisp);
   userDefinedVel  = (double *) dbg_alloca(sizeof(double)*claw->numUserDisp);
   userSupFunc->usd_disp( midtime, userDefinedDisp, userDefinedVel );
   int i;
   for(i=0; i<claw->numUserDisp; ++i) {
  
     int dsaDof = dsa->locate(claw->userDisp[i].nnum, 
                         1 << claw->userDisp[i].dofnum);
     if(dsaDof < 0)
       continue;
     int cdof = c_dsa->invRCN(dsaDof);
     if(cdof >= 0)
      state.prescDisp[cdof] = userDefinedDisp[i];
   }
 
   // KENDALL: pass coordinate set to update Prescribed displacements
   // this way the control file can be uniformized with the linear one!
 
   // domain->getNodes());
   setBC(userDefinedDisp, userDefinedVel);
 }
 
 elementInternalForce.zero();
 getStiffAndForce(state, residual, elementInternalForce, *totalRes);
                       //   kelArray, residual);
 times->buildStiffAndForce +=  getTime();
 
 // return residual force norm
 return residual.norm();
}

void
NLMatProbDesc::formRHSinitializer(Vector &fext, Vector &velocity, Vector &elementInternalForce,
                                  NLState &st, Vector &rhs, NLState *refState)
{
  // rhs = (fext - fint - Cv)
  rhs = fext;
  // PJSA fint and Cv assumed zero for now
}

void
NLMatProbDesc::formRHSpredictor(Vector &velocity, Vector &residual, Vector &rhs,
   NLState &st, double) 
{
  times->predictorTime -= getTime();
  velocity.print("Vel");
  residual.print("FExt");
 
  if(temp == 0)
    temp = new Vector(rhs);
 
  Vector &temp1 = *temp;
 
  // Vector temp( rhs );
  temp1.zero();
 
  temp1.linC( st.disp, delta, velocity );
  // rhs = M*(disp + delta*velocity)
  M->mult( temp1, rhs );
 
  // rhs = M*(disp + delta*velocity) + delta^2*residual
  rhs.linC( 1.0, rhs, delta*delta, residual );
 
  times->predictorTime += getTime();
}

void
NLMatProbDesc::addInertialTerm(Vector &glRes, Vector &acc)
{
  Muc->mult(*boundAcc, *unconstrInertiaF);
  Mcc->mult(*boundAcc, *boundInertiaF); 
  Muc->transposeMultAdd(acc.data(),  boundInertiaF->data());
  M->multAdd(acc.data(), unconstrInertiaF->data());
  
  int *uMap = c_dsa->getUnconstrNum();
  int *cMap = c_dsa->getConstrndNum();
  for(int i = 0; i < dsa->size(); ++i) {
    if(uMap[i] >= 0)
      glRes[i] -= (*unconstrInertiaF)[uMap[i]];
    else
      glRes[i] -= (*boundInertiaF)[cMap[i]];
  }
  
}
int
NLMatProbDesc::checkConvergence(int iteration, double normRes,
                                Vector& residual, Vector& dv,
                                double time)
{
     /*if(dofTypeArray == 0)
       dofTypeArray = cdsa->makeDofTypeArray();*/
#ifdef PRINT_FORCENORMS 
     ConstrainedDSA *cdsa = domain.getCDSA();
     double momenNorm = 0.0;
     double forceNorm = 0.0;
     int i;
     for(i=0; i<domain.numNodes(); ++i) {
       int dof = cdsa->locate(i, DofSet::Xdisp);
       if(dof >= 0)
         forceNorm += (residual[dof]*residual[dof]);
       dof = cdsa->locate(i, DofSet::Ydisp);
       if(dof >= 0)
         forceNorm += (residual[dof]*residual[dof]);
       dof = cdsa->locate(i, DofSet::Zdisp);
       if(dof >= 0)
         forceNorm += (residual[dof]*residual[dof]);
       dof = cdsa->locate(i, DofSet::Xrot);
       if(dof >= 0)
         momenNorm += (residual[dof]*residual[dof]);
       dof = cdsa->locate(i, DofSet::Yrot);
       if(dof >= 0)
         momenNorm += (residual[dof]*residual[dof]);
       dof = cdsa->locate(i, DofSet::Zrot);
       if(dof >= 0)
         momenNorm += (residual[dof]*residual[dof]);
     }
 
     if(iteration == 0)  {
       firstForceNorm = sqrt(forceNorm);
       firstMomenNorm = sqrt(momenNorm);
       if(firstForceNorm == 0) firstForceNorm = 1;
       if(firstMomenNorm == 0) firstMomenNorm = 1;
     }
     //if(verboseFlag)
       fprintf(stderr,"===> time %f force Norm %e Moment Norm %e\n",
               time,
               sqrt(forceNorm)/firstForceNorm,
               sqrt(momenNorm)/firstMomenNorm);
#endif 

     double normDv     = dv.norm();
     double normEnergy = residual*dv;
     if(iteration == 0)  {
       firstRes = normRes;
       firstDv  = normDv;
       firstEng = normEnergy;
     }

     double relRes = (firstRes == 0) ? 0 : normRes/firstRes;
     double relDv  = (firstDv == 0) ? 0 : normDv /firstDv;
     double relEng = (firstEng == 0) ? 0 : normEnergy/firstEng;
 
     int converged = 0;
 
     if( relRes <= tolerance)
       converged = 1;
 
     // Check for divergence
     if(relRes >= 10000 || relDv >= 10000)
       converged = -1;
 
     if(verboseFlag) {
       fprintf(stderr," Iteration # %d\n",iteration);
       fprintf(stderr," r      = %e dv      = %e energy      = %e\n"
                      " rel. r = %e rel. dv = %e rel. energy = %e\n",
                        normRes,normDv,normEnergy,
                        relRes,relDv,relEng);
     }
 
     totIter++;
// XXX XML
//     fprintf(res,"%d %e %e %e %e %e\n",totIter,time,normRes,relRes, normDv, relDv);
//     fflush(res);
 
     // Store residual norm and dv norm for output
// XXX XML
/*
     times->norms[numSystems].normDv      = normDv;
     times->norms[numSystems].relativeDv  = relDv;
     times->norms[numSystems].normRes     = normRes;
     times->norms[numSystems].relativeRes = relRes;
     times->numSystems = numSystems;
*/

     numSystems += 1;
 
     return converged;
}

void
NLMatProbDesc::reBuild(NLState &, int iteration, double delta)
{
 times->rebuild -= getTime();

 // Rebuild every updateK iterations
 if( iteration % domain.solInfo().getNLInfo().updateK == 0 ) {
/* PJSA
   solver->reBuild(kelArray, melArray, delta);
*/
   spm->zeroAll();
   AllOps<double> ops;
   ops.Kuc = kuc;
   double beta, gamma, alphaf, alpham, dt = 2*delta;
   getNewmarkParameters(beta, gamma, alphaf, alpham);
   double Kcoef = dt*dt*beta;
   double Ccoef = dt*gamma;
   double Mcoef = (1-alpham)/(1-alphaf);
   domain.makeSparseOps<double>(ops, Kcoef, Mcoef, Ccoef, spm, kelArray, melArray);
   solver->factor();
 }

 times->rebuild += getTime();

}

int
NLMatProbDesc::reBuild(int iteration, int, NLState &)
{
 times->rebuild -= getTime();

 //cerr << "here in NLMatProbDesc::reBuild(int iteration, int, NLState &)\n";
 // Rebuild every updateK iterations
 if( iteration % domain.solInfo().getNLInfo().updateK == 0 ) {
/* PJSA
   solver->reBuild(kelArray);
*/
   spm->zeroAll();
   AllOps<double> ops;
   domain.makeSparseOps<double>(ops, 1.0, 0.0, 0.0, spm, kelArray, (FullSquareMatrix *) NULL);
   solver->factor();
 }

 times->rebuild += getTime();

 return (iteration % domain.solInfo().getNLInfo().updateK == 0) ? 1 : 0;
}

double
NLMatProbDesc::formRHScorrector(Vector &inc_displacement, Vector &velocity,
                                Vector &residual,         Vector &rhs)
{
 times->correctorTime -= getTime();
 
 if(temp == 0)
   temp = new Vector(rhs);
 
 Vector &temp1 = *temp;
 
 // Vector temp( rhs );
 temp1.zero();
 
 temp1.linC( inc_displacement, -delta, velocity );
 
 // rhs = M*temp1
 M->mult( temp1, rhs );
 
 // rhs = delta^2*residual - M(inc_displacement - delta*velocity)
 rhs.linC( delta*delta, residual, -1.0, rhs );
 
 times->correctorTime += getTime();   
 return residual.norm(); 
}

void
NLMatProbDesc::printTimers(double& timeLoop)  
{
 //long memoryUsed = solver->size();
 //double solveTime     = solver->getSolutionTime();
 
 //if(solver->numRBM() != 0)
 //   fprintf(stderr," ... Number of rbm %4d             ...\n",
 //           solver->numRBM());
 
 // times->printStaticTimers( solveTime, memoryUsed, domain, timeLoop );
}

void
NLMatProbDesc::setBC(double *, double *)
{
 // MLX We need to look at what this function should do.
}

void
NLMatProbDesc::extract(Vector &, Vector &, Vector &, double *, double *, double *)
{
 // MLX We need to look at what this function should do.
}

void
NLState::get_inc_displacement (Vector &inc, NLState &prev, bool)
{
 inc = disp;
 inc -= prev.disp;
}

void
NLState::midpoint_step_update(Vector &vel_n, double delta, NLState &ss, Vector &acc_n)
{
 int size = disp.size();
 double coef = 2.0/delta;
 
 int i;
 for(i = 0; i < size; ++i)
   vel_n[i] = coef*(disp[i] - ss.disp[i])-vel_n[i];
 for(i = 0; i < size; ++i) {
   double t = disp[i];
   disp[i] = 2*disp[i]-ss.disp[i];
   acc_n[i] = 0.25*coef*coef*(disp[i]-ss.disp[i])-coef*vel_n[i]-acc_n[i];
   ss.disp[i] = t;
 }

 // Now we need to update the internal state of the material
}

void
NLState::printNode(int)
{

}

/*void
NLState::update(Vector &)
{

}

void
NLState::get_inc_displacement(Vector &, NLState &)
{

}

void
NLState::midpoint_step_update(Vector &, double &, NLState &)
{

}*/

void
NLMatProbDesc::getNewmarkParameters(double &beta, double &gamma,
                                    double &alphaf, double &alpham)
{
 beta  = domain.solInfo().newmarkBeta;
 gamma = domain.solInfo().newmarkGamma;
 alphaf = domain.solInfo().newmarkAlphaF;
 alpham = domain.solInfo().newmarkAlphaM;
}

