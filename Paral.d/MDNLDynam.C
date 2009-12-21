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
  // AERO not supported
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

  return domain->solInfo().aeroFlag;
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
  // rhs = delta^2 * residual - M * (inc_displacement - delta * velocity) - C * delta * inc_displacement
  localTemp->linC(inc_displacement, -localDelta, velocity);
  M->mult(*localTemp, rhs);
  rhs.linC(localDelta * localDelta, residual, -1.0, rhs);
  if(C) { // DAMPING
    C->mult(inc_displacement, *localTemp);
    rhs.linAdd(-localDelta, *localTemp);
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

  // rhs = M*velocity
  M->mult(velocity, rhs);

  // rhs = delta*M*velocity + delta^2*residual
  rhs.linC(localDelta, rhs, localDelta * localDelta, residual);

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
                              DistrVector& elementInternalForce, double midtime) 
{
 // note: midtime is used with claw (currently not supported)
 times->buildStiffAndForce -= getTime();

 if(midtime != -1.0 && claw && userSupFunc) {
   if(claw->numUserDisp > 0) {
     double *userDefineDisp = new double[claw->numUserDisp];
     double *userDefineVel  = new double[claw->numUserDisp];
     userSupFunc->usd_disp(midtime, userDefineDisp, userDefineVel);
     execParal2R(decDomain->getNumSub(), this, &MDNLDynamic::subUpdateGeomStateUSDD, geomState, userDefineDisp);
     delete [] userDefineDisp; delete [] userDefineVel;
   }
 }

 execParal3R(decDomain->getNumSub(), this, &MDNLDynamic::subGetStiffAndForce, geomState,
             residual, elementInternalForce);

 updateMpcRhs(geomState);

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

 if(iteration % domain->solInfo().getNLInfo().updateK == 0) { // XXXX
   //filePrint(stderr,"===> REBUILDING TANGENT STIFFNESS MATRIX\n");
   times->norms[numSystems].rebuildTang = 1;
   execParal(decDomain->getNumSub(), this, &MDNLDynamic::rebuildKelArray, localDelta);
   solver->reBuild(kelArray, geomState, iteration); 
 } else
   times->norms[numSystems].rebuildTang = 0;

 times->rebuild += getTime();
}

void
MDNLDynamic::rebuildKelArray(int isub, double localDelta)
{
  SubDomain *sd = decDomain->getSubDomain(isub);
  Connectivity *allDofs = sd->getAllDOFs();
  double delta2 = localDelta * localDelta;
  SparseMatrix *kuc = (Kuc) ? (*Kuc)[isub] : 0;
  if(kuc) kuc->zeroAll();
  for(int iele = 0; iele < sd->numElements(); ++iele) {
    int dim = kelArray[isub][iele].dim();
    if(kuc) kuc->add(kelArray[isub][iele], (*allDofs)[iele]);
    for(int i = 0; i < dim; ++i)
      for(int j = 0; j < dim; ++j) {
        if(C) kelArray[isub][iele][i][j] = delta2 * kelArray[isub][iele][i][j] + localDelta * celArray[isub][iele][i][j] + melArray[isub][iele][i][j];
        else kelArray[isub][iele][i][j] = delta2 * kelArray[isub][iele][i][j] + melArray[isub][iele][i][j];
      }
  }
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
  decDomain->buildOps(*allOps, Kcoef, Mcoef, Ccoef);
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
MDNLDynamic::getExternalForce(DistrVector& f, DistrVector& gravityForce,
                              int tIndex, double time, DistrGeomState* geomState,
                              DistrVector& elementInternalForce, DistrVector &aeroF)
{
  // AERO and Thermal currently not supported
  times->formRhs -= getTime();

  // add the FORCE (including MFTT), HDNB, ROBIN, GRAVITY and PRESSURE forces
  // XXXX make sure geomState is updated at time t (for follower pressure) including USDD
  execParal5R(decDomain->getNumSub(), this, &MDNLDynamic::subGetExternalForce,
              f, gravityForce, *geomState, tIndex, time);

  // add the USDF forces
  if(claw && userSupFunc) {
    if(claw->numUserForce > 0) {
      double *userDefineForce = new double[claw->numUserForce];
      userSupFunc->usd_forc(time, userDefineForce);
      decDomain->addUserForce(f, userDefineForce);
      delete [] userDefineForce; 
    }
  }

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
      for(int i=0; i<decDomain->getNumSub(); ++i) subExtractControlDisp(i, *geomState, ctrdisp);
#ifdef DISTRIBUTED
      structCom->globalMax(claw->numSensor, ctrdisp);
#endif

      userSupFunc->ctrl(ctrdisp, ctrvel, ctracc, ctrfrc, time);
      decDomain->addCtrl(f, ctrfrc);

      delete [] ctrdisp; delete [] ctrvel; delete [] ctracc; delete [] ctrfrc;
    }
  }

/* aero is not supported yet for multidomain
  // get solver info
  SolverInfo& sinfo = domain->solInfo();
  if(sinfo.aeroFlag >= 0) {

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
    if(aero_f) {
      aero_f->zero();
      aero_f->linAdd(alpha, *aeroForce, (1.0-alpha), *prevFrc);
    }

    *prevFrc= *aeroForce;
    prevTime = tFluid;
    prevIndex = tIndex;
  }
*/
  times->formRhs += getTime();
}

void
MDNLDynamic::subGetExternalForce(int isub, DistrVector& f, DistrVector& gravityForce, 
                                 DistrGeomState& geomState, int tIndex, double time)
{
  StackVector localf(f.subData(isub),f.subLen(isub));
  StackVector localg(gravityForce.subData(isub), gravityForce.subLen(isub));

  SubDomain *sd = decDomain->getSubDomain(isub);

  PrevFrc dummy(0); // used for aero 
  sd->computeExtForce4(dummy, localf, localg, tIndex, time);

  // add the PRESSURE forces
  if(sd->pressureFlag()) sd->buildPressureForce(localf, geomState[isub]);
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
   if(aeroAlg >= 0) cerr << "MDNLDynamic::dynamCommToFluid is not implemented\n";
}

