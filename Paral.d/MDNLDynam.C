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

  MultiDomainOp mdop(&MultiDomainOp::getConstForce, decDomain->getAllSubDomains(), &v);
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
  sd->updatePrescribedDisp(geomState[isub], getDeltaLambda());
}

void
MDNLDynamic::formRHSinitializer(DistrVector &fext, DistrVector &velocity, DistrVector &elementInternalForce, 
                                  DistrGeomState &geomState, DistrVector &rhs)
{
  // rhs = (fext - fint - Cv)
  rhs = fext;
  elementInternalForce.zero();
  getStiffAndForce(geomState, rhs, elementInternalForce);
  if(C) {
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
  addMpcForces(rhs); // XXXX

  double resN = sqrt(solver->getFNormSq(rhs));
  times->correctorTime += getTime();
  //return rhs.norm();
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

  updateMpcRhs(geomState); // XXXX

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
  dt        = domain->solInfo().dt;
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

  execParal3R(decDomain->getNumSub(), this, &MDNLDynamic::subGetStiffAndForce, geomState,
              residual, elementInternalForce);

  updateMpcRhs(geomState);

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
                                 DistrVector &res, DistrVector &elemIntForce)
{
  // PJSA: 10-4-2007 copied from MDNLStatic
  SubDomain *sd = decDomain->getSubDomain(isub);
  StackVector residual(res.subData(isub), res.subLen(isub));
  // eIF = element internal force
  StackVector eIF(elemIntForce.subData(isub), elemIntForce.subLen(isub));
  sd->getStiffAndForce(*geomState[isub], eIF, allCorot[isub], kelArray[isub], residual);

  //sd->updateMpcRhs(*geomState[isub]);
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
   double Mcoef = (1-alpham)/(1-alphaf);
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
  solver = (FetiSolver *) allOps->dynMat;
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

  if(domain->solInfo().aeroFlag >= 0)
    cerr << " *** WARNING: aeroelastic is not supported for Multi-domain nonlinear dynamics \n";

  localTemp = new DistrVector(decDomain->solVecInfo());

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
  if(sinfo.aeroFlag >= 0) {

    double gamma = sinfo.newmarkGamma;
    double alphaf = sinfo.newmarkAlphaF;

    aeroForce->zero();
    int iscollocated;
    double tFluid = distFlExchanger->getFluidLoad(*aeroForce, tIndex, t,
                                                  alphaf, iscollocated);
    if (iscollocated == 0) {
      if(prevIndex >= 0) {
        *aeroForce *= (1/gamma);
        aeroForce->linAdd(((gamma-1.0)/gamma),*prevFrc);
      }
    }

    double alpha = 1.0-alphaf;
    if(prevIndex < 0) alpha = 1.0;

    f.linAdd(alpha, *aeroForce, (1.0-alpha), *prevFrc);
    aero_f.zero();
    aero_f.linAdd(alpha, *aeroForce, (1.0-alpha), *prevFrc);

    *prevFrc = *aeroForce;
    prevTime = tFluid;
    prevIndex = tIndex;
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

FetiSolver *
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
#ifdef DISTRIBUTED
  // DistDom output doesn't support veloc & acc 
  decDomain->postProcessing(geomState, allCorot, time);
#else
  SysState<DistrVector> distState(ext_force, vel_n, acc_n, vel_p); 
  decDomain->postProcessing(geomState, allCorot, time, &distState);
#endif
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

void
MDNLDynamic::initNewton()
{
  solver->initNewton();
}

void
MDNLDynamic::addMpcForces(DistrVector& vec)
{
  execParal1R(decDomain->getNumSub(), this, &MDNLDynamic::subAddMpcForces, vec);
}

void
MDNLDynamic::subAddMpcForces(int isub, DistrVector& vec)
{
  SubDomain *sd = decDomain->getSubDomain(isub);
  double *mpcForces = new double[sd->numMPCs()];      // don't delete  
  solver->getLocalMpcForces(isub, mpcForces);         // mpcForces set to incremental mpc lagrange multipliers
  sd->addMpcForceIncrement(mpcForces);                // mpcForces set to total mpc lagrange multipliers
  StackVector localvec(vec.subData(isub), vec.subLen(isub));
  sd->constraintProductTmp(mpcForces, localvec);      // C^T*lambda added to vec
}

void
MDNLDynamic::updateMpcRhs(DistrGeomState &geomState)
{
  execParal1R(decDomain->getNumSub(), this, &MDNLDynamic::subUpdateMpcRhs, geomState);
}

void
MDNLDynamic::subUpdateMpcRhs(int isub, DistrGeomState &geomState)
{
  SubDomain *sd = decDomain->getSubDomain(isub);
  sd->updateMpcRhs(*geomState[isub], decDomain->getMpcToSub());
}

void 
MDNLDynamic::dynamCommToFluid(DistrGeomState* geomState, DistrGeomState* bkGeomState,
                              DistrVector& velocity, DistrVector& bkVelocity,
                              DistrVector& vp, DistrVector& bkVp, int step, int parity,
                              int aeroAlg) 
{  
  if(domain->solInfo().aeroFlag >= 0 && !domain->solInfo().lastIt) {
    DistrVector d_n_aero(decDomain->solVecInfo()); d_n_aero.zero();
    execParal5R(decDomain->getNumSub(), this, &MDNLDynamic::subDynamCommToFluid, d_n_aero, geomState, bkGeomState, parity, aeroAlg);
    if(!parity && aeroAlg == 5) {
      velocity.linC(0.5, velocity, 0.5, bkVelocity);
      vp.linC(0.5, vp, 0.5, bkVp);
    }
    DistrVector acceleration(decDomain->solVecInfo()); acceleration.zero(); // XXXX

    SysState<DistrVector> state(d_n_aero, velocity, acceleration, vp);

    if (verboseFlag)
      filePrint(stderr," ... Send displacements to Fluid at step %d\n",(step+1));
    distFlExchanger->sendDisplacements(state, usrDefDisps, usrDefVels);
  }
}

void
MDNLDynamic::subDynamCommToFluid(int isub, DistrVector& v, DistrGeomState* distrGeomState, 
                                 DistrGeomState* bkDistrGeomState, int parity, int aeroAlg)
{
  SubDomain *sd = decDomain->getSubDomain(isub);
  StackVector d_n(v.subData(isub), v.subLen(isub));
  Vector d_n2(v.subLen(isub), 0.0);
  GeomState* geomState = (*distrGeomState)[isub];
  GeomState* bkGeomState = (*bkDistrGeomState)[isub];
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
    int numDofs = dsa[iSub]->numNodes();
    usrDefDisps[iSub] = new double[numDofs];
    usrDefVels[iSub] = new double[numDofs];

    for(int iDof = 0; iDof < numDofs; iDof++)  {
      usrDefDisps[iSub][iDof] = 0.0;
      usrDefVels[iSub][iDof] = 0.0;
    }
  }

  // create distributed fluid exchanger
  distFlExchanger = new DistFlExchanger(cs, elemSet, cdsa, dsa);
  //mddPostPro->setPostProcessor(distFlExchanger);
  //mddPostPro->setUserDefs(usrDefDisps, usrDefVels);

  // negotiate with the fluid code
  distFlExchanger->negotiate();

  int restartinc = (sinfo.nRestart >= 0) ? (sinfo.nRestart) : 0;

  //fprintf(stderr,"struct cpu %d Sending Parameters to Fluid Code\n", structCom->myID());

  distFlExchanger->sendParam(sinfo.aeroFlag, sinfo.dt, sinfo.tmax, restartinc,
                             sinfo.isCollocated, sinfo.alphas);

  // initialize the Parity
  if(sinfo.aeroFlag == 5 || sinfo.aeroFlag == 4) {
     distFlExchanger->initRcvParity(1);
     distFlExchanger->initSndParity(1);
   } else {
     distFlExchanger->initRcvParity(-1);
     distFlExchanger->initSndParity(-1);
   }

  // create distributed state vector
  SysState<DistrVector> state(disp, vel, accel, lastVel);

  // send initial displacements
  distFlExchanger->sendDisplacements(state, usrDefDisps, usrDefVels);

  return sinfo.aeroFlag;
}

int
MDNLDynamic::getAeroAlg()
{
  return domain->solInfo().aeroFlag;
}

void
MDNLDynamic::thermoePreProcess()
{
  filePrint(stderr, "Paral.d/MDNLDynam.C: thermoePreProcess not implemented here\n");
}

int
MDNLDynamic::getThermoeFlag()
{
  return domain->solInfo().thermoeFlag;
}
