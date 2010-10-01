#include <stdio.h>

#include <Threads.d/Paral.h>
#include <Driver.d/Domain.h>
#include <Paral.d/MDNLDynam.h>
#include <Math.d/SparseMatrix.h>
#include <Solvers.d/Solver.h>
#include <Driver.d/SubDomain.h>
#include <Threads.d/PHelper.h>
#include <Timers.d/StaticTimers.h>
#include <Timers.d/GetTime.h>
#include <Math.d/Vector.h>
#include <Math.d/mathUtility.h>
#include <Utils.d/DistHelper.h>
#include <Paral.d/MDStatic.h>
#include <Paral.d/MDDynam.h>
#ifdef DISTRIBUTED
#include <Dist.d/DistDom.h>
#endif
#include <Hetero.d/DistFlExchange.h>
#include <Utils.d/ModeData.h>

extern ModeData modeData;

// ***************************************************************
// *                                                             *
// *  Purpose: Multiple Domain implementation of nonlinear       *
// *           dynamics				                 *	
// *                                                             *
// *  Coded by: Kendall H. Pierson and Philip J. S. Avery        *
// ***************************************************************

void
MDNLDynamic::getConstForce(DistrVector &v)
{
  times->formRhs -= getTime();

  MultiDomainOp mdop(&MultiDomainOp::getConstForce, decDomain->getAllSubDomains(), &v, Kuc);
  threadManager->execParal(decDomain->getNumSub(), &mdop);

  times->formRhs += getTime();
}

int
MDNLDynamic::getInitState(DistrVector &d_n, DistrVector &v_n, DistrVector &a_n, DistrVector &v_p)
{
  MultiDomainOp mdop(&MultiDomainOp::getInitState, decDomain->getAllSubDomains(),
                     &d_n, &v_n, &a_n, &v_p);
  threadManager->execParal(decDomain->getNumSub(), &mdop);

  if(claw && userSupFunc && claw->numSensor) {
    double *ctrdisp = new double[claw->numSensor];
    double *ctrvel  = new double[claw->numSensor];
    double *ctracc  = new double[claw->numSensor];
#ifdef DISTRIBUTED
    for(int i=0; i<claw->numSensor; ++i) ctrdisp[i] = ctrvel[i] = ctracc[i] = std::numeric_limits<double>::min();
#endif
    decDomain->extractControlData(d_n, v_n, a_n, ctrdisp, ctrvel, ctracc);
#ifdef USE_MPI
    structCom->globalMax(claw->numSensor, ctrdisp);
    structCom->globalMax(claw->numSensor, ctrvel);
    structCom->globalMax(claw->numSensor, ctracc);
#endif
    userSupFunc->init(ctrdisp, ctrvel, ctracc);
    delete [] ctrdisp; delete [] ctrvel; delete [] ctracc;
  }

  int aeroAlg = domain->solInfo().aeroFlag;
  // call aeroPreProcess if a restart file does not exist
  if(aeroAlg >= 0 && geoSource->getCheckFileInfo()->lastRestartFile == 0)
    aeroPreProcess(d_n, v_n, a_n, v_p);

  if(domain->solInfo().thermoeFlag >= 0)
    thermoePreProcess();

  if(domain->solInfo().thermohFlag >= 0)
    thermohPreProcess(d_n);

  if(domain->solInfo().aeroheatFlag >= 0)
    aeroheatPreProcess(d_n, v_n, v_p);

  return aeroAlg;
}

void
MDNLDynamic::updatePrescribedDisplacement(DistrGeomState *geomState)
{
  times->timePresc -= getTime();
  execParal1R(decDomain->getNumSub(), this, &MDNLDynamic::subUpdatePrescribedDisplacement, *geomState);
  times->timePresc += getTime();
}

void
MDNLDynamic::subUpdatePrescribedDisplacement(int isub, DistrGeomState& geomState)
{
  SubDomain *sd = decDomain->getSubDomain(isub);
  sd->updatePrescribedDisp(geomState[isub]);
}

void
MDNLDynamic::formRHSinitializer(DistrVector &fext, DistrVector &velocity, DistrVector &elementInternalForce, 
                                  DistrGeomState &geomState, DistrVector &rhs)
{
  // rhs = (fext - fint - Cv)
  rhs = fext;
  elementInternalForce.zero();
  getStiffAndForce(geomState, rhs, elementInternalForce);
  if(domain->solInfo().order == 2 && C) {
    C->mult(velocity, *localTemp);
    rhs.linC(rhs, -1.0, *localTemp);
  }
}

double
MDNLDynamic::formRHScorrector(DistrVector& inc_displacement, DistrVector& velocity, DistrVector& acceleration,
                              DistrVector& residual, DistrVector& rhs, double localDelta)
{
  // PJSA 10-4-2007 copied from NonLinDynamic
  times->correctorTime -= getTime();
  if(domain->solInfo().order == 1) {
    M->mult(inc_displacement, rhs);
    rhs.linC(localDelta, residual, -1.0, rhs);
  }
  else {
    double beta, gamma, alphaf, alpham, dt = 2*localDelta;
    getNewmarkParameters(beta, gamma, alphaf, alpham);
    // rhs = dt*dt*beta*residual - ((1-alpham)/(1-alphaf)*M+dt*gamma*C)*inc_displacement
    //       + (dt*(1-alpham)*M - dt*dt*(beta-(1-alphaf)*gamma)*C)*velocity
    //       + (dt*dt*((1-alpham)/2-beta)*M - dt*dt*dt*(1-alphaf)*(2*beta-gamma)/2*C)*acceleration
    localTemp->linC(-(1-alpham)/(1-alphaf), inc_displacement, dt*(1-alpham), velocity, dt*dt*((1-alpham)/2-beta), acceleration);
    M->mult(*localTemp, rhs);
    if(C) {
      localTemp->linC(-dt*gamma, inc_displacement, -dt*dt*(beta-(1-alphaf)*gamma), velocity, -dt*dt*dt*(1-alphaf)*(2*beta-gamma)/2, acceleration);
      C->multAdd(*localTemp, rhs);
    }
    rhs.linAdd(dt*dt*beta, residual);
  }

  double resN = sqrt(solver->getFNormSq(rhs));
  times->correctorTime += getTime();
  return resN;
}

void
MDNLDynamic::formRHSpredictor(DistrVector& velocity, DistrVector& acceleration, DistrVector& residual,
                              DistrVector& rhs, DistrGeomState &geomState,
                              double midtime, double localDelta)
{
  // PJSA 10-4-2007 copied from single domain equivalent in Problems.d/NonLinDynamic.C
  times->predictorTime -= getTime();

  if(claw && userSupFunc) {
    if(claw->numUserDisp > 0) {

      // allocate memory for user defined motion
      double *userDefineDisp = new double[claw->numUserDisp];
      double *userDefineDispLast = new double[claw->numUserDisp];
      double *userDefineVel = new double[claw->numUserDisp];

      // get user defined motion
      userSupFunc->usd_disp(midtime, userDefineDisp, userDefineVel);
      userSupFunc->usd_disp(midtime-delta, userDefineDispLast, userDefineVel);

      // update state
      execParal2R(decDomain->getNumSub(), this, &MDNLDynamic::subUpdateGeomStateUSDD, geomState, userDefineDisp);

      // get delta disps
      for(int j = 0; j < claw->numUserDisp; j++)
        userDefineDisp[j] -= userDefineDispLast[j];

      // update force residual with KUC
      if(Kuc) 
        execParal2R(decDomain->getNumSub(), this, &MDNLDynamic::subKucTransposeMultSubtractClaw, residual, userDefineDisp);

      delete [] userDefineDisp; delete [] userDefineDispLast; delete [] userDefineVel;
    }
  }

  if(domain->solInfo().order == 1)
    rhs.linC(residual, localDelta);
  else {
    double beta, gamma, alphaf, alpham, dt = 2*localDelta;
    getNewmarkParameters(beta, gamma, alphaf, alpham);
    // rhs = dt*dt*beta*residual + (dt*(1-alpham)*M - dt*dt*(beta-(1-alphaf)*gamma)*C)*velocity
    //       + (dt*dt*((1-alpham)/2-beta)*M - dt*dt*dt*(1-alphaf)*(2*beta-gamma)/2*C)*acceleration
    localTemp->linC(dt*(1-alpham), velocity, dt*dt*((1-alpham)/2-beta), acceleration);
    M->mult(*localTemp, rhs);
    if(C) {
      localTemp->linC(-dt*dt*(beta-(1-alphaf)*gamma), velocity, -dt*dt*dt*(1-alphaf)*(2*beta-gamma)/2, acceleration);
      C->multAdd(*localTemp, rhs);
    }
    rhs.linAdd(dt*dt*beta, residual);
  }

  times->predictorTime += getTime();
}

void
MDNLDynamic::subKucTransposeMultSubtractClaw(int isub, DistrVector& residual, double *userDefineDisp)
{
  SubDomain *sd = decDomain->getSubDomain(isub);
  ControlLawInfo *subClaw = sd->getClaw();
  if(subClaw) {
    if(subClaw->numUserDisp) {
      double *subUserDefineDisp = new double[subClaw->numUserDisp];
      for(int i=0; i<subClaw->numUserDisp; ++i) {
        subUserDefineDisp[i] = userDefineDisp[sd->getUserDispDataMap()[i]];
      }
      (*Kuc)[isub]->transposeMultSubtractClaw(subUserDefineDisp, residual.subData(isub), subClaw->numUserDisp, clawDofs[isub]);
      delete [] subUserDefineDisp;
    }
  }
}

void
MDNLDynamic::computeTimeInfo()
{
  // Time integration information

  // Get total time and time step size and store them 
  totalTime = domain->solInfo().tmax;
  dt        = domain->solInfo().getTimeStep();
  delta     = 0.5*dt;

  // Compute maximum number of steps
  maxStep = (int) ( (totalTime+0.49*dt)/dt );

  // Compute time remainder
  double remainder = totalTime - maxStep*dt;
  if(std::abs(remainder)>0.01*dt){
    domain->solInfo().tmax = totalTime = maxStep*dt;
    fprintf(stderr, " Warning: Total time is being changed to : %e\n", totalTime);
  }

  // XXXX
  // set half time step size in user defined functions 
  //if(userSupFunc)
  //  userSupFunc->setDt(delta);

}

MDNLDynamic::MDNLDynamic(Domain *d)
{ 
  domain = d;
#ifdef DISTRIBUTED
  decDomain = new GenDistrDomain<double>(domain);
#else
  decDomain = new GenDecDomain<double>(domain);
#endif
  numSystems = 0;
  secondRes = 0.0;
  claw = 0; 
  userSupFunc = 0;
  mu = 0; lambda = 0;
}

MDNLDynamic::~MDNLDynamic()
{
  if(mu) delete [] mu;
  if(lambda) delete [] lambda;
}

int
MDNLDynamic::checkConvergence(int iteration, double normRes, DistrVector &residual, DistrVector& dv, double time)
{
  times->timeCheck -= getTime();
  double normDv  = dv.norm();
  double normEnergy = residual*dv;

  if(iteration == 0)  {
    firstRes = normRes;
    firstDv  = normDv;
    firstEng = normEnergy;
  }

  if(iteration == 1)
    secondRes = normRes;

  double relRes = normRes/firstRes;
  double relDv  = normDv /firstDv;

  int converged = 0;
  // Check for convergence
  if(normRes <= tolerance*firstRes) converged = 1;
  // Check for divergence
  else if(normRes >= 1.0e10 * firstRes && normRes > secondRes) converged = -1;

#ifdef PRINT_RESIDUALS
  double relEng = normEnergy/firstEng;
  if(verboseFlag) {
    fprintf(stderr," Iteration # %d\n",iteration);
    fprintf(stderr," r      = %e dv      = %e energy      = %e\n"
                   " rel. r = %e rel. dv = %e rel. energy = %e\n",
                     normRes,normDv,normEnergy,
                     relRes,relDv,relEng);
  }
#endif
  totIter++;

  // Store residual norm and dv norm for output
  times->norms[numSystems].normDv      = normDv;
  times->norms[numSystems].relativeDv  = relDv;
  times->norms[numSystems].normRes     = normRes;
  times->norms[numSystems].relativeRes = relRes;
  times->numSystems = numSystems;

  numSystems += 1;
  times->timeCheck += getTime();
  return converged;
}

double
MDNLDynamic::getStiffAndForce(DistrGeomState& geomState, DistrVector& residual,
                              DistrVector& elementInternalForce, double t) 
{
  times->buildStiffAndForce -= getTime();

  // update the geomState according to the USDD prescribed displacements
  if(claw && userSupFunc) {
    if(claw->numUserDisp > 0) {
      double *userDefineDisp = new double[claw->numUserDisp];
      double *userDefineVel  = new double[claw->numUserDisp];
      userSupFunc->usd_disp(t, userDefineDisp, userDefineVel);
      execParal2R(decDomain->getNumSub(), this, &MDNLDynamic::subUpdateGeomStateUSDD, geomState, userDefineDisp);
      delete [] userDefineDisp; delete [] userDefineVel;
    }
  }

  execParal4R(decDomain->getNumSub(), this, &MDNLDynamic::subGetStiffAndForce, geomState,
              residual, elementInternalForce, t);

  if(t != -1.0) updateConstraintTerms(&geomState);

  // add the ACTUATOR forces
  if(claw && userSupFunc) {
    if(claw->numActuator > 0) {
      double *ctrdisp = new double[claw->numSensor];
      double *ctrvel  = new double[claw->numSensor];
      double *ctracc  = new double[claw->numSensor];
      double *ctrfrc  = new double[claw->numActuator];

      for(int i=0; i<claw->numSensor; ++i) ctrvel[i] = ctracc[i] = 0.0; // not supported
#ifdef DISTRIBUTED
      for(int i=0; i<claw->numSensor; ++i) ctrdisp[i] = std::numeric_limits<double>::min();
#endif
      for(int i=0; i<decDomain->getNumSub(); ++i) subExtractControlDisp(i, geomState, ctrdisp);
#ifdef DISTRIBUTED
      structCom->globalMax(claw->numSensor, ctrdisp);
#endif

      userSupFunc->ctrl(ctrdisp, ctrvel, ctracc, ctrfrc, t);
      decDomain->addCtrl(residual, ctrfrc);

      delete [] ctrdisp; delete [] ctrvel; delete [] ctracc; delete [] ctrfrc;
    }
  }

  times->buildStiffAndForce += getTime();

  return sqrt(solver->getFNormSq(residual));
}

void
MDNLDynamic::subGetStiffAndForce(int isub, DistrGeomState &geomState,
                                 DistrVector &res, DistrVector &elemIntForce, double t)
{
  // PJSA: 10-4-2007 copied from MDNLStatic
  SubDomain *sd = decDomain->getSubDomain(isub);
  StackVector residual(res.subData(isub), res.subLen(isub));
  // eIF = element internal force
  StackVector eIF(elemIntForce.subData(isub), elemIntForce.subLen(isub));
  sd->getStiffAndForce(*geomState[isub], eIF, allCorot[isub], kelArray[isub], residual, 1.0, t);
}

void
MDNLDynamic::subUpdateGeomStateUSDD(int isub, DistrGeomState &geomState, double *userDefineDisp)
{
  SubDomain *sd = decDomain->getSubDomain(isub);
  ControlLawInfo *subClaw = sd->getClaw();
  if(subClaw) {
    if(subClaw->numUserDisp) {
      double *subUserDefineDisp = new double[subClaw->numUserDisp];
      for(int i=0; i<subClaw->numUserDisp; ++i)
        subUserDefineDisp[i] = userDefineDisp[sd->getUserDispDataMap()[i]];
      geomState[isub]->updatePrescribedDisplacement(subUserDefineDisp, subClaw, sd->getNodes());
      delete [] subUserDefineDisp;
    }
  }
}

DistrGeomState*
MDNLDynamic::createGeomState()
{
  times->timeGeom -= getTime();
  return new DistrGeomState(decDomain);
  times->timeGeom += getTime();
}

DistrGeomState*
MDNLDynamic::copyGeomState(DistrGeomState* geomState)
{
 times->timeGeom -= getTime();
 return new DistrGeomState(*geomState);
 times->timeGeom += getTime();
}

void
MDNLDynamic::reBuild(DistrGeomState& geomState, int iteration, double localDelta)
{
 times->rebuild -= getTime();

 if(iteration % domain->solInfo().getNLInfo().updateK == 0) {
   times->norms[numSystems].rebuildTang = 1;
   GenMDDynamMat<double> ops;
   ops.sysSolver = solver;
   ops.Kuc = Kuc;
   double beta, gamma, alphaf, alpham, dt = 2*localDelta;
   getNewmarkParameters(beta, gamma, alphaf, alpham);
   double Kcoef = (domain->solInfo().order == 1) ? localDelta : dt*dt*beta;
   double Ccoef = (domain->solInfo().order == 1) ? 0.0 : dt*gamma;
   double Mcoef = (domain->solInfo().order == 1) ? 1 : (1-alpham)/(1-alphaf);
   decDomain->rebuildOps(ops, Mcoef, Ccoef, Kcoef, kelArray, melArray);
 } else
   times->norms[numSystems].rebuildTang = 0;

 times->rebuild += getTime();
}

void
MDNLDynamic::preProcess()
{
  // Structure used to store timers
  times = new StaticTimers;

  times->memoryPreProcess -= threadManager->memoryUsed();

  totIter = 0;
  
  // Set the nonlinear tolerance
  tolerance = domain->solInfo().getNLInfo().tolRes;

  // Constructs all renumbering, connectivities and dofsets
  times->preProcess -= getTime();
  decDomain->preProcess();
  times->preProcess += getTime();

  // Make each subdomain's dofs
  times->makeDOFs -= getTime();
  execParal(decDomain->getNumSub(), this, &MDNLDynamic::makeSubDofs);
  times->makeDOFs += getTime();

  // Make each subdomain's corotators
  times->corotatorTime -= getTime();
  allCorot = new Corotator**[decDomain->getNumSub()]; 
  execParal(decDomain->getNumSub(), this, &MDNLDynamic::makeSubCorotators);
  times->corotatorTime += getTime();

  times->memoryPreProcess += threadManager->memoryUsed();

  // Construct FETI Solver
  times->getFetiSolverTime -= getTime();
  allOps = new MDDynamMat;
  double Kcoef = 0.0;
  double Mcoef = 1.0;
  double Ccoef = 0.0;
  decDomain->buildOps(*allOps, Mcoef, Ccoef, Kcoef);
  solver = (ParallelSolver *) allOps->dynMat;
  M = allOps->M;
  C = (domain->solInfo().alphaDamp != 0 || domain->solInfo().betaDamp != 0) ? allOps->C : 0;
  Kuc = allOps->Kuc;
  times->getFetiSolverTime += getTime();

  times->memoryPreProcess -= threadManager->memoryUsed();

  // Allocate each subdomain's array of stiffness matrices
  kelArray = new FullSquareMatrix*[decDomain->getNumSub()];

  // Allocate each subdomain's array of mass matrices
  melArray = new FullSquareMatrix*[decDomain->getNumSub()];

  if(C) celArray = new FullSquareMatrix*[decDomain->getNumSub()];

  // Now make those arrays
  execParal(decDomain->getNumSub(), this, &MDNLDynamic::makeSubElementArrays);

  // Look if there is a user supplied routine for control
  claw = geoSource->getControlLaw();

  // create list of usdd node dofs mapped to cdsa dof numbers
  if(claw && claw->numUserDisp)  {
    clawDofs = new int * [decDomain->getNumSub()];
    execParal(decDomain->getNumSub(), this, &MDNLDynamic::makeSubClawDofs);
  }
  // Check to see if there is a user supplied function
  // for displacements, forces or control law
  userSupFunc = domain->getUserSuppliedFunction();

  localTemp = new DistrVector(decDomain->solVecInfo());

  domain->InitializeStaticContactSearch(MortarHandler::CTC, decDomain->getNumSub(), decDomain->getAllSubDomains());
  mu = new std::map<int,double>[decDomain->getNumSub()];
  lambda = new std::vector<double>[decDomain->getNumSub()];

  times->memoryPreProcess -= threadManager->memoryUsed();
}

void
MDNLDynamic::makeSubDofs(int isub)
{
  SubDomain *sd = decDomain->getSubDomain(isub);
  sd->makeAllDOFs();
}

void
MDNLDynamic::makeSubCorotators(int isub)
{
  SubDomain *sd  = decDomain->getSubDomain(isub);
  int numele     = sd->numElements();
  allCorot[isub] = new Corotator*[numele];
  sd->createCorotators(allCorot[isub]);
}

void
MDNLDynamic::makeSubElementArrays(int isub)
{
  SubDomain *sd = decDomain->getSubDomain(isub);
  if(C) sd->createKelArray(kelArray[isub], melArray[isub], celArray[isub]); 
  else sd->createKelArray(kelArray[isub], melArray[isub]);
}

void
MDNLDynamic::makeSubClawDofs(int isub)
{
  SubDomain *sd = decDomain->getSubDomain(isub);
  ControlLawInfo *subClaw = sd->getClaw();
  int nClaw = subClaw->numUserDisp;
  if(nClaw > 0) {
    clawDofs[isub] = new int[nClaw];
    for(int j = 0; j < nClaw; ++j) {
      int dd = sd->getDSA()->locate(subClaw->userDisp[j].nnum, (1 << subClaw->userDisp[j].dofnum));
      clawDofs[isub][j] = sd->getCDSA()->invRCN(dd);
    }
  }
  else clawDofs[isub] = 0;
}

void
MDNLDynamic::getExternalForce(DistrVector& f, DistrVector& constantForce,
                              int tIndex, double t, DistrGeomState* geomState,
                              DistrVector& elementInternalForce, DistrVector& aero_f)
{
  times->formRhs -= getTime();

  // update nodal temperature for thermoe
  if(domain->solInfo().thermoeFlag >= 0 && tIndex >= 0) {
    distFlExchanger->getStrucTemp(nodalTemps->data());
    if(verboseFlag) filePrint(stderr," ... [E] Received temperatures ...\n");
  }

  execParal3R(decDomain->getNumSub(), this, &MDNLDynamic::subGetExternalForce,
              f, constantForce, t);

  // add the USDF forces
  if(claw && userSupFunc) {
    if(claw->numUserForce > 0) {
      double *userDefineForce = new double[claw->numUserForce];
      userSupFunc->usd_forc(t, userDefineForce);
      decDomain->addUserForce(f, userDefineForce);
      delete [] userDefineForce; 
    }
  }

  // AERO
  SolverInfo& sinfo = domain->solInfo();
  if(sinfo.aeroFlag >= 0 && tIndex >= 0) {

    double gamma = sinfo.newmarkGamma;
    double alphaf = sinfo.newmarkAlphaF;

    aeroForce->zero();
    int iscollocated;
    double tFluid = distFlExchanger->getFluidLoad(*aeroForce, tIndex, t,
                                                  alphaf, iscollocated);
    if(verboseFlag) filePrint(stderr," ... [E] Received fluid forces ...\n");
    if (iscollocated == 0) {
      if(prevIndex >= 0) {
        *aeroForce *= (1/gamma);
        aeroForce->linAdd(((gamma-1.0)/gamma),*prevFrc);
      }
    }

    double alpha = 1.0-alphaf;
    if(prevIndex < 0) alpha = 1.0;

    aero_f.linC(alpha, *aeroForce, (1.0-alpha), *prevFrc);
    f += aero_f;

    *prevFrc = *aeroForce;
    prevTime = tFluid;
    prevIndex = tIndex;
  }

  // AEROH
  if(sinfo.aeroheatFlag >= 0 && tIndex >= 0) {

    aeroForce->zero();
    double tFluid = distFlExchanger->getFluidFlux(*aeroForce, tIndex, t);
    if(verboseFlag) filePrint(stderr," ... [T] Received fluid fluxes (%e) ...\n", aeroForce->sqNorm());

    /*  Compute fluid flux at n+1/2, since we use midpoint rule in thermal */

    int useProjector = domain->solInfo().filterFlags;

    if(tIndex == 0)
      f += *aeroForce;
    else {
      if(useProjector) f = *aeroForce;
      else
        f.linAdd(0.5, *aeroForce, 0.5, *prevFrc);
    }

    *prevFrc = *aeroForce;
  }

  times->formRhs += getTime();
}

void
MDNLDynamic::subGetExternalForce(int isub, DistrVector& f, DistrVector& constantForce, double time)
{
  StackVector localf(f.subData(isub), f.subLen(isub));
  StackVector localg(constantForce.subData(isub), constantForce.subLen(isub));

  SubDomain *sd = decDomain->getSubDomain(isub);

  sd->computeExtForce4(localf, localg, time);
}

void
MDNLDynamic::subExtractControlDisp(int isub, DistrGeomState &geomState, double *ctrdsp)
{
  SubDomain *sd = decDomain->getSubDomain(isub);
  ControlLawInfo *subClaw = sd->getClaw();
  if(subClaw && subClaw->numSensor > 0) {
    CoordSet &nodes = sd->getNodes();
    NodeState *nodeState = geomState[isub]->getNodeState();
    for(int i=0; i<subClaw->numSensor; ++i) {
      switch (subClaw->sensor[i].dofnum) {
        case 0:
          ctrdsp[sd->getSensorDataMap()[i]] = nodeState[subClaw->sensor[i].nnum].x
                        - nodes[subClaw->sensor[i].nnum]->x;
          break;
        case 1:
          ctrdsp[sd->getSensorDataMap()[i]] = nodeState[subClaw->sensor[i].nnum].y
                        - nodes[subClaw->sensor[i].nnum]->y;
          break;
        case 2:
          ctrdsp[sd->getSensorDataMap()[i]] = nodeState[subClaw->sensor[i].nnum].z
                        - nodes[subClaw->sensor[i].nnum]->z;
          break;
        default:
          fprintf(stderr, "ERROR: Sensor dof %d not available in MDNLDynam::subExtractControlDisp\n",subClaw->sensor[i].dofnum+1);
      }
    }
  }
}

ParallelSolver *
MDNLDynamic::getSolver()
{
  return solver;
}

MultiDomainPostProcessor *
MDNLDynamic::getPostProcessor()
{
  return new MultiDomainPostProcessor(decDomain, solver);
}

void
MDNLDynamic::printTimers(double timeLoop)
{
  int i;
  long (*memory)=(long *) dbg_alloca(sizeof(long)*decDomain->getNumSub());
  for (i = 0; i < decDomain->getNumSub(); ++i)
    memory[i] = 0;

  MultiDomainPostProcessor *mdpp = getPostProcessor();

  execParal(decDomain->getNumSub(), mdpp,
           &MultiDomainPostProcessor::getMemoryPrec, memory);
  long totMemPrec = 0;
  for(i = 0; i < decDomain->getNumSub(); ++i)
    totMemPrec += memory[i];

  for(i = 0; i < decDomain->getNumSub(); ++i)
    memory[i] = 0;

  execParal(decDomain->getNumSub(), mdpp,
           &MultiDomainPostProcessor::getMemoryK, memory);

  long totMemK = 0;
  for(i=0; i<decDomain->getNumSub(); ++i)
    totMemK += memory[i];

  Timings &fetiTimers = solver->getTimers();
  fetiTimers.preconditioner.addOverAll(totMemPrec, 0.0);
  fetiTimers.kMem.addOverAll(totMemK, 0.0);

#ifdef DISTRIBUTED

  double mem1 = (double) totMemPrec;
  if(structCom) mem1 = structCom->globalSum(mem1);
  totMemPrec = (long) mem1;

  mem1 = (double) totMemK;
  if(structCom) mem1 = structCom->globalSum(mem1);
  totMemK = (long) mem1;

#endif

  times->memoryK = totMemK;
  times->memoryPrecond = totMemPrec;

  times->timeTimers -= getTime();

  times->numSubdomain = decDomain->getNumSub();

  times->printTimers(domain, solver->getTimers(),
                     solver->getSolutionTime());

  times->timeTimers += getTime();
}

void
MDNLDynamic::dynamOutput(DistrGeomState *geomState, DistrVector &vel_n, DistrVector &vel_p, 
                         double time, int index, DistrVector &ext_force, DistrVector &aeroF, DistrVector &acc_n)
{
  if(claw && claw->numUserDisp) {
    double *userDefineDisp = new double[claw->numUserDisp];
    double *userDefineVel  = new double[claw->numUserDisp];
    userSupFunc->usd_disp(time,userDefineDisp,userDefineVel);
    execParal2R(decDomain->getNumSub(), this, &MDNLDynamic::subUpdateGeomStateUSDD, *geomState, userDefineDisp);
    paralApply(decDomain->getNumSub(), decDomain->getAllSubDomains(), &GenSubDomain<double>::setUserDefBC, userDefineDisp, userDefineVel);
    delete [] userDefineDisp; delete [] userDefineVel;
  }
  SysState<DistrVector> distState(ext_force, vel_n, acc_n, vel_p); 
  decDomain->postProcessing(geomState, allCorot, time, &distState);
}

void
MDNLDynamic::processLastOutput()  
{
  OutputInfo *oinfo = geoSource->getOutputInfo();
  for(int iOut = 0; iOut < geoSource->getNumOutInfo(); iOut++)
    oinfo[iOut].interval = 1;
}

void
MDNLDynamic::getInitialTime(int &initTimeIndex, double &initTime)
{
  initTimeIndex = domain->solInfo().initialTimeIndex;
  initTime      = domain->solInfo().initialTime;
}

void
MDNLDynamic::getNewmarkParameters(double &beta, double &gamma,
                                  double &alphaf, double &alpham)
{
  beta  = domain->solInfo().newmarkBeta;
  gamma = domain->solInfo().newmarkGamma;
  alphaf = domain->solInfo().newmarkAlphaF;
  alpham = domain->solInfo().newmarkAlphaM;
}

int
MDNLDynamic::getMaxit()
{
  return domain->solInfo().getNLInfo().maxiter;
}

double
MDNLDynamic::getDeltaLambda()
{
  return domain->solInfo().getNLInfo().dlambda;
}

int
MDNLDynamic::getNumStages()
{
  return int(0.2 + domain->solInfo().getNLInfo().maxLambda / domain->solInfo().getNLInfo().dlambda);
}

DistrInfo&
MDNLDynamic::solVecInfo()
{
  return decDomain->solVecInfo();
}

DistrInfo&
MDNLDynamic::elemVecInfo()
{ 
  return *decDomain->elementVectorInfo();
} 
  
DistrInfo&
MDNLDynamic::sysVecInfo()
{ 
  return decDomain->sysVecInfo();
}

double
MDNLDynamic::getResidualNorm(DistrVector &r)
{
 //returns: sqrt( (r+c^T*lambda)**2 + pos_part(gap)**2 )
 DistrVector w(r);
 execParal1R(decDomain->getNumSub(), this, &MDNLDynamic::addConstraintForces, w); // w = r + C^T*lambda
                  // note C = grad(gap) has already been updated in getStiffAndForce.
                  // XXXX need to make sure lambda_i is correctly mapped to C_i. I think this is done
                  // correctly only for the case of one contactsurfaces pair
 return sqrt(solver->getFNormSq(w));
}

void
MDNLDynamic::addConstraintForces(int isub, DistrVector& vec)
{
  // I need to treat the contact forces from CONTACTSURFACES separately due to search,
  // the ith lagrange multiplier at iteration n may not correspond to the ith constraint
  // after updating the contact surfaces
  SubDomain *sd = decDomain->getSubDomain(isub);
  StackVector localvec(vec.subData(isub), vec.subLen(isub));
  sd->addConstraintForces(mu[isub], lambda[isub], localvec);  // C^T*lambda added to vec
}

void
MDNLDynamic::getConstraintMultipliers(int isub)
{
  SubDomain *sd = decDomain->getSubDomain(isub);
  mu[isub].clear();
  lambda[isub].clear();
  sd->getConstraintMultipliers(mu[isub], lambda[isub]);
}

void
MDNLDynamic::updateConstraintTerms(DistrGeomState* geomState)
{
  GenFetiDPSolver<double> *fetiSolver = dynamic_cast<GenFetiDPSolver<double> *>(solver);
  if(fetiSolver) {
    execParal(decDomain->getNumSub(), this, &MDNLDynamic::getConstraintMultipliers);
    if(domain->GetnContactSurfacePairs()) {
      // this function updates the linearized contact conditions (the lmpc coeffs are the gradient and the rhs is the gap)
      // XXXX the hessian of the constraint functions needs to be computed also
      domain->UpdateSurfaces(MortarHandler::CTC, geomState, decDomain->getAllSubDomains());
      domain->PerformStaticContactSearch(MortarHandler::CTC); // note: dynamic contact search not supported by acme for face-face interactions
      domain->deleteSomeLMPCs(mpc::ContactSurfaces);
      domain->ExpComputeMortarLMPC(MortarHandler::CTC);
      domain->CreateMortarToMPC();
      decDomain->reProcessMPCs();
      fetiSolver->reconstructMPCs(decDomain->mpcToSub_dual, decDomain->mpcToMpc, decDomain->mpcToCpu);
    }
    // set the gap for the linear constraints
    decDomain->setConstraintGap(geomState, fetiSolver);
  }
}

void 
MDNLDynamic::dynamCommToFluid(DistrGeomState* geomState, DistrGeomState* bkGeomState,
                              DistrVector& velocity, DistrVector& bkVelocity,
                              DistrVector& vp, DistrVector& bkVp, int step, int parity,
                              int aeroAlg) 
{  
  if(domain->solInfo().aeroFlag >= 0 && !domain->solInfo().lastIt) {
    DistrVector d_n(decDomain->solVecInfo()); d_n.zero();
    execParal5R(decDomain->getNumSub(), this, &MDNLDynamic::subDynamCommToFluid, d_n, geomState, bkGeomState, parity, aeroAlg);
    if(!parity && aeroAlg == 5) {
      velocity.linC(0.5, velocity, 0.5, bkVelocity);
      vp.linC(0.5, vp, 0.5, bkVp);
    }
    DistrVector acceleration(decDomain->solVecInfo()); acceleration.zero(); // XXXX

    SysState<DistrVector> state(d_n, velocity, acceleration, vp);

    distFlExchanger->sendDisplacements(state, usrDefDisps, usrDefVels);
    if(verboseFlag) filePrint(stderr," ... [E] Sent displacements to Fluid at step %d ...\n", step+1);
  }

  if(domain->solInfo().aeroheatFlag >= 0) {
    DistrVector d_n(decDomain->solVecInfo()); d_n.zero();
    execParal2R(decDomain->getNumSub(), this, &MDNLDynamic::subDynamCommToFluidAeroheat, d_n, geomState);

    SysState<DistrVector> tempState(d_n, velocity, vp);

    distFlExchanger->sendTemperature(tempState);
    if(verboseFlag) filePrint(stderr," ... [T] Sent temperatures to Structure at step %d ...\n", step+1);
  }

  if(domain->solInfo().thermohFlag >= 0) {
    for(int i = 0; i < decDomain->getNumSub(); ++i) {
      SubDomain *sd = decDomain->getSubDomain(i);
      for(int j = 0; j < sd->numNodes(); ++j) {
        nodalTemps->subData(i)[j] = (*(*geomState)[i])[j].x;
      }
    }

    distFlExchanger->sendStrucTemp(*nodalTemps);
    if(verboseFlag) filePrint(stderr," ... [T] Sent temperatures to Structure at step %d ...\n", step+1);
  }
}

void
MDNLDynamic::subDynamCommToFluid(int isub, DistrVector& v, DistrGeomState* distrGeomState, 
                                 DistrGeomState* bkDistrGeomState, int parity, int aeroAlg)
{
  SubDomain *sd = decDomain->getSubDomain(isub);
  StackVector d_n(v.subData(isub), v.subLen(isub));
  GeomState* geomState = (*distrGeomState)[isub];
  double* bcx = usrDefDisps[isub];

  // Make d_n_aero from geomState
  ConstrainedDSA *c_dsa = sd->getCDSA();
  DofSetArray *dsa = sd->getDSA();
  CoordSet &nodes = sd->getNodes();
  int numNodes = sd->numNodes();

  for(int i = 0; i < numNodes; ++i) {

    int xloc  = c_dsa->locate(i, DofSet::Xdisp );
    int xloc1 =   dsa->locate(i, DofSet::Xdisp );

    if(xloc >= 0)
      d_n[xloc]  = ( (*geomState)[i].x - nodes[i]->x);
    else if (xloc1 >= 0)
      bcx[xloc1] = ( (*geomState)[i].x - nodes[i]->x);

    int yloc  = c_dsa->locate(i, DofSet::Ydisp );
    int yloc1 =   dsa->locate(i, DofSet::Ydisp );

    if(yloc >= 0)
      d_n[yloc]  = ( (*geomState)[i].y - nodes[i]->y);
    else if (yloc1 >= 0)
      bcx[yloc1] = ( (*geomState)[i].y - nodes[i]->y);

    int zloc  = c_dsa->locate(i, DofSet::Zdisp);
    int zloc1 =   dsa->locate(i, DofSet::Zdisp);

    if(zloc >= 0)
      d_n[zloc]  = ( (*geomState)[i].z - nodes[i]->z);
    else if (zloc1 >= 0)
      bcx[zloc1] = ( (*geomState)[i].z - nodes[i]->z);
  }

  if(!parity && aeroAlg == 5) {
    Vector d_n2(v.subLen(isub), 0.0);
    GeomState* bkGeomState = (*bkDistrGeomState)[isub];
    for(int i = 0; i < numNodes; ++i) {

      int xloc  = c_dsa->locate(i, DofSet::Xdisp );
      if(xloc >= 0)
        d_n2[xloc]  = ( (*bkGeomState)[i].x - nodes[i]->x);

      int yloc  = c_dsa->locate(i, DofSet::Ydisp );
      if(yloc >= 0)
        d_n2[yloc]  = ( (*bkGeomState)[i].y - nodes[i]->y);

      int zloc  = c_dsa->locate(i, DofSet::Zdisp);
      if(zloc >= 0)
        d_n2[zloc]  = ( (*bkGeomState)[i].z - nodes[i]->z);
    }
    d_n.linC(0.5, d_n, 0.5, d_n2);
  }
}

void
MDNLDynamic::subDynamCommToFluidAeroheat(int isub, DistrVector& v, DistrGeomState* distrGeomState)
{
  SubDomain *sd = decDomain->getSubDomain(isub);
  StackVector d_n(v.subData(isub), v.subLen(isub));
  GeomState* geomState = (*distrGeomState)[isub];
  double* bcx = usrDefDisps[isub];

  // Make d_n from geomState
  ConstrainedDSA *c_dsa = sd->getCDSA();
  DofSetArray *dsa = sd->getDSA();
  CoordSet &nodes = sd->getNodes();
  int numNodes = sd->numNodes();

  for(int i = 0; i < numNodes; ++i) {

    int xloc  = c_dsa->locate(i, DofSet::Temp );
    int xloc1 =   dsa->locate(i, DofSet::Temp );

    if(xloc >= 0)
      d_n[xloc]  = (*geomState)[i].x;
    else if (xloc1 >= 0)
      bcx[xloc1] = (*geomState)[i].x;
  }
}

void 
MDNLDynamic::readRestartFile(DistrVector &d_n, DistrVector &v_n, DistrVector &a_n,
                             DistrVector &v_p, DistrGeomState &geomState)
{
  ControlInfo *cinfo = geoSource->getCheckFileInfo();
  if(cinfo->lastRestartFile)
    filePrint(stderr, "Paral.d/MDNLDynam.C: readRestartFile not implemented here\n");
}

int
MDNLDynamic::aeroPreProcess(DistrVector &disp, DistrVector &vel,
                            DistrVector &accel, DistrVector &lastVel)
{
  // get solver info
  SolverInfo& sinfo = domain->solInfo();

  // Initialize previous force data
  // Reexamine for the case of restart
  prevFrc = new DistrVector(solVecInfo());
  prevIndex = -1;
  prevTime = 0;

  // Initialize the aeroforce vector
  aeroForce = new DistrVector(solVecInfo());

  if(sinfo.aeroFlag < 0)
    return 0;

  Connectivity *cpuToSub = geoSource->getCpuToSub();

  // get cpu id
#ifdef USE_MPI
  int myId = structCom->myID();
#else
  int myId = 0;
#endif

  int numLocSub = cpuToSub->num(myId);

  SubDomain **subdomain = decDomain->getAllSubDomains();

  // allocate for pointer arrays
  CoordSet **cs = new CoordSet *[numLocSub];
  Elemset **elemSet = new Elemset *[numLocSub];
  DofSetArray **cdsa = new DofSetArray *[numLocSub];
  DofSetArray **dsa = new DofSetArray *[numLocSub];
  usrDefDisps = new double *[numLocSub];
  usrDefVels = new double *[numLocSub];

  int iSub;
  for(iSub = 0; iSub < numLocSub; iSub++)  {

    // assemble coordsets in this mpi
    cs[iSub] = &subdomain[iSub]->getNodes();

    // assemble element sets in this mpi
    elemSet[iSub] = &subdomain[iSub]->getElementSet();

    // assemble constrained and unconstrained dofset arrays in this mpi
    cdsa[iSub] = subdomain[iSub]->getCDSA();
    dsa[iSub] = subdomain[iSub]->getDSA();

    // allocate and initialize for the user defined disps and vels
    int numDofs = dsa[iSub]->size();
    usrDefDisps[iSub] = new double[numDofs];
    usrDefVels[iSub] = new double[numDofs];

    for(int iDof = 0; iDof < numDofs; iDof++)  {
      usrDefDisps[iSub][iDof] = 0.0;
      usrDefVels[iSub][iDof] = 0.0;
    }
  }

  int numOutInfo = geoSource->getNumOutInfo();
  OutputInfo *oinfo = geoSource->getOutputInfo();

  int flag = 0;

  // Check if aero forces are requested for output
  int iInfo;
  for(iInfo = 0; iInfo < numOutInfo; ++iInfo) {
    if(oinfo[iInfo].type == OutputInfo::AeroForce) {
      flag = 1;
      break;
    }
  }

  // create distributed fluid exchanger
  if(flag)
    distFlExchanger = new DistFlExchanger(cs, elemSet, cdsa, dsa, oinfo+iInfo);
  else
    distFlExchanger = new DistFlExchanger(cs, elemSet, cdsa, dsa);

  // negotiate with the fluid code
  distFlExchanger->negotiate();

  int restartinc = (sinfo.nRestart >= 0) ? (sinfo.nRestart) : 0;

  DistrVector dispAero(disp);

  if(sinfo.gepsFlg == 1) {
    // If we are in the first time step, and we initialized with
    // IDISP6, do not send IDISP6
    if(domain->numInitDisp() == 0 && sinfo.zeroInitialDisp != 1) {
      filePrint(stderr," ... DO NOT SEND IDISP6 0\n"); //HB
    } else {
      filePrint(stderr," ... SENDING IDISP6 0\n"); //HB
      for(iSub = 0; iSub < numLocSub; iSub++) {
        BCond* iDis6 = subdomain[iSub]->getInitDisp6();
        for(int i = 0; i < subdomain[iSub]->numInitDisp6(); ++i) {
          int dof = cdsa[iSub]->locate(iDis6[i].nnum, 1 << iDis6[i].dofnum);
          if(dof >= 0)
            dispAero[dof] += iDis6[i].val;
        }
      }
    }
  }

  SysState<DistrVector> state(dispAero, vel, accel, lastVel);

  if(sinfo.aeroFlag == 8) {
    distFlExchanger->sendParam(sinfo.aeroFlag, sinfo.getTimeStep(), sinfo.mppFactor,
                               restartinc, sinfo.isCollocated, sinfo.alphas);
    distFlExchanger->sendModeFreq(modeData.frequencies, modeData.numModes);
    if(verboseFlag) filePrint(stderr,"... [E] Sent parameters and mode frequencies ...\n");
    distFlExchanger->sendModeShapes(modeData.numModes, modeData.numNodes,
                                    modeData.modes, state, sinfo.mppFactor);
    if(verboseFlag) filePrint(stderr,"... [E] Sent mode shapes ...\n");
  }
  else {
    distFlExchanger->sendParam(sinfo.aeroFlag, sinfo.getTimeStep(), sinfo.tmax, restartinc,
                               sinfo.isCollocated, sinfo.alphas);
    if(verboseFlag) filePrint(stderr,"... [E] Sent parameters ...\n");

    // initialize the Parity
    if(sinfo.aeroFlag == 5 || sinfo.aeroFlag == 4) {
       distFlExchanger->initRcvParity(1);
       distFlExchanger->initSndParity(1);
     } else {
       distFlExchanger->initRcvParity(-1);
       distFlExchanger->initSndParity(-1);
     }

    // send initial displacements
    distFlExchanger->sendDisplacements(state, usrDefDisps, usrDefVels);
    if(verboseFlag) filePrint(stderr,"... [E] Sent initial displacements ...\n");

    if(sinfo.aeroFlag == 1) { // Ping pong only
      filePrint(stderr, "Ping Pong Only requested. Structure code exiting\n");
    }
  }

  return sinfo.aeroFlag;
}

void
MDNLDynamic::thermoePreProcess()
{
  if(domain->solInfo().thermoeFlag >=0) {

    Connectivity *cpuToSub = geoSource->getCpuToSub();
    int myId = structCom->myID();
    int numLocSub = cpuToSub->num(myId);
    SubDomain **subdomain = decDomain->getAllSubDomains();

    // if sinfo.aeroFlag >= 0, flExchanger has already been initialize before,
    // thus, only when sinfo.aeroFlag < 0 is necessary.
    if(domain->solInfo().aeroFlag < 0) {

      // allocate for pointer arrays
      CoordSet **cs = new CoordSet *[numLocSub];
      Elemset **elemSet = new Elemset *[numLocSub];
      DofSetArray **cdsa = new DofSetArray *[numLocSub];
      DofSetArray **dsa = new DofSetArray *[numLocSub];
      usrDefDisps = new double *[numLocSub];
      usrDefVels = new double *[numLocSub];

      int iSub;
      for(iSub = 0; iSub < numLocSub; iSub++)  {

        // assemble coordsets in this mpi
        cs[iSub] = &subdomain[iSub]->getNodes();

        // assemble element sets in this mpi
        elemSet[iSub] = &subdomain[iSub]->getElementSet();

        // assemble constrained and unconstrained dofset arrays in this mpi
        cdsa[iSub] = subdomain[iSub]->getCDSA();
        dsa[iSub] = subdomain[iSub]->getDSA();

        // allocate and initialize for the user defined disps and vels
        int numDofs = dsa[iSub]->size();
        usrDefDisps[iSub] = new double[numDofs];
        usrDefVels[iSub] = new double[numDofs];

        for(int iDof = 0; iDof < numDofs; iDof++)  {
          usrDefDisps[iSub][iDof] = 0.0;
          usrDefVels[iSub][iDof] = 0.0;
        }
      }

      distFlExchanger = new DistFlExchanger(cs, elemSet, cdsa, dsa);
      //mddPostPro->setPostProcessor(distFlExchanger);
      //mddPostPro->setUserDefs(usrDefDisps, usrDefVels);
    }

    nodalTemps = new DistrVector(decDomain->ndVecInfo());
    for(int iSub = 0; iSub < numLocSub; iSub++) subdomain[iSub]->temprcvd = nodalTemps->subData(iSub); // XXXX
    int buffLen = nodalTemps->size();

    distFlExchanger->thermoread(buffLen);

    distFlExchanger->getStrucTemp(nodalTemps->data()) ;
    if(verboseFlag) filePrint(stderr," ... [E] Received initial temperatures ...\n");
  }
}

void
MDNLDynamic::thermohPreProcess(DistrVector& d)
{
  if(domain->solInfo().thermohFlag >=0) {

    Connectivity *cpuToSub = geoSource->getCpuToSub();
    int myId = structCom->myID();
    int numLocSub = cpuToSub->num(myId);
    SubDomain **subdomain = decDomain->getAllSubDomains();

    // if sinfo.aeroheatFlag >= 0, flExchanger has already been initialize before,
    // thus, only when sinfo.aeroheatFlag < 0 is necessary.
    if(domain->solInfo().aeroheatFlag < 0) {

      // allocate for pointer arrays
      CoordSet **cs = new CoordSet *[numLocSub];
      Elemset **elemSet = new Elemset *[numLocSub];
      DofSetArray **cdsa = new DofSetArray *[numLocSub];
      DofSetArray **dsa = new DofSetArray *[numLocSub];
      usrDefDisps = new double *[numLocSub];
      usrDefVels = new double *[numLocSub];

      int iSub;
      for(iSub = 0; iSub < numLocSub; iSub++)  {

        // assemble coordsets in this mpi
        cs[iSub] = &subdomain[iSub]->getNodes();

        // assemble element sets in this mpi
        elemSet[iSub] = &subdomain[iSub]->getElementSet();

        // assemble constrained and unconstrained dofset arrays in this mpi
        cdsa[iSub] = subdomain[iSub]->getCDSA();
        dsa[iSub] = subdomain[iSub]->getDSA();

        // allocate and initialize for the user defined disps and vels
        int numDofs = dsa[iSub]->size();
        usrDefDisps[iSub] = new double[numDofs];
        usrDefVels[iSub] = new double[numDofs];

        for(int iDof = 0; iDof < numDofs; iDof++)  {
          usrDefDisps[iSub][iDof] = 0.0;
          usrDefVels[iSub][iDof] = 0.0;
        }
      }

      distFlExchanger = new DistFlExchanger(cs, elemSet, cdsa, dsa);
      //mddPostPro->setPostProcessor(distFlExchanger);
      //mddPostPro->setUserDefs(usrDefDisps, usrDefVels);
    }

    nodalTemps = new DistrVector(decDomain->ndVecInfo());
    int buffLen = nodalTemps->size();
    //mddPostPro->setNodalTemps(nodalTemps);

    distFlExchanger->thermoread(buffLen);

    for(int i = 0; i < decDomain->getNumSub(); ++i) {
      SubDomain *sd = decDomain->getSubDomain(i);
      for(int j = 0; j < sd->numNodes(); ++j) {
        int tloc  = sd->getCDSA()->locate(j, DofSet::Temp);
        int tloc1 = sd->getDSA()->locate(j, DofSet::Temp);
        double temp  = (tloc >= 0) ? d.subData(i)[tloc] : sd->getBcx()[tloc1];
        if(tloc1 < 0) temp = 0.0;
        nodalTemps->subData(i)[j] = temp;
      }
    }

    distFlExchanger->sendStrucTemp(*nodalTemps);
    if(verboseFlag) filePrint(stderr," ... [T] Sent initial temperatures ...\n");
  }

}

void
MDNLDynamic::aeroheatPreProcess(DistrVector &disp, DistrVector &vel, DistrVector &lastVel)
{
  // get solver info
  SolverInfo& sinfo = domain->solInfo();

  // Initialize previous force data
  // Reexamine for the case of restart
  prevFrc = new DistrVector(solVecInfo());
  prevIndex = -1;
  prevTime = 0;

  // Initialize the aeroforce vector
  aeroForce = new DistrVector(solVecInfo());

  if(sinfo.aeroheatFlag < 0)
    return;

  Connectivity *cpuToSub = geoSource->getCpuToSub();

  // get cpu id
#ifdef USE_MPI
  int myId = structCom->myID();
#else
  int myId = 0;
#endif

  int numLocSub = cpuToSub->num(myId);

  SubDomain **subdomain = decDomain->getAllSubDomains();

  // allocate for pointer arrays
  CoordSet **cs = new CoordSet *[numLocSub];
  Elemset **elemSet = new Elemset *[numLocSub];
  DofSetArray **cdsa = new DofSetArray *[numLocSub];
  DofSetArray **dsa = new DofSetArray *[numLocSub];
  usrDefDisps = new double *[numLocSub];
  usrDefVels = new double *[numLocSub];

  int iSub;
  for(iSub = 0; iSub < numLocSub; iSub++)  {

    // assemble coordsets in this mpi
    cs[iSub] = &subdomain[iSub]->getNodes();

    // assemble element sets in this mpi
    elemSet[iSub] = &subdomain[iSub]->getElementSet();

    // assemble constrained and unconstrained dofset arrays in this mpi
    cdsa[iSub] = subdomain[iSub]->getCDSA();
    dsa[iSub] = subdomain[iSub]->getDSA();

    // allocate and initialize for the user defined disps
    int numDofs = dsa[iSub]->size();
    usrDefDisps[iSub] = new double[numDofs];

    for(int iDof = 0; iDof < numDofs; iDof++)  {
      usrDefDisps[iSub][iDof] = 0.0;
    }
  }

  int numOutInfo = geoSource->getNumOutInfo();
  OutputInfo *oinfo = geoSource->getOutputInfo();

  int flag = 0;

  // Check if aero fluxes are requested for output
  int iInfo;
  for(iInfo = 0; iInfo < numOutInfo; ++iInfo) {
    if(oinfo[iInfo].type == OutputInfo::AeroForce) {
      flag = 1;
      break;
    }
  }

  // create distributed fluid exchanger
  if(flag)
    distFlExchanger = new DistFlExchanger(cs, elemSet, cdsa, dsa, oinfo+iInfo);
  else
    distFlExchanger = new DistFlExchanger(cs, elemSet, cdsa, dsa);

  // negotiate with the fluid code
  distFlExchanger->negotiate();

  int restartinc = (sinfo.nRestart >= 0) ? (sinfo.nRestart) : 0;

  SysState<DistrVector> state(disp, vel, lastVel);

  distFlExchanger->sendTempParam(sinfo.aeroheatFlag, sinfo.getTimeStep(), sinfo.tmax, restartinc,
                                 sinfo.alphat);
  if(verboseFlag) filePrint(stderr,"... [T] Sent parameters ...\n");

  // send initial displacements
  distFlExchanger->sendTemperature(state);
  if(verboseFlag) filePrint(stderr,"... [T] Sent initial temperatures ...\n");
}

int
MDNLDynamic::getAeroAlg()
{
  return domain->solInfo().aeroFlag;
}

int
MDNLDynamic::getThermoeFlag()
{
  return domain->solInfo().thermoeFlag;
}

int
MDNLDynamic::getThermohFlag()
{
  return domain->solInfo().thermohFlag;
}

int
MDNLDynamic::getAeroheatFlag()
{
  return domain->solInfo().aeroheatFlag;
}
