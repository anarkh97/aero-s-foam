#include <iostream>
#include <Driver.d/Domain.h>
#include <Driver.d/DynamProbType.h>
#include <Paral.d/MDDynam.h>
#include <Driver.d/Dynam.h>
#include <Paral.d/MDOp.h>
#include <Timers.d/StaticTimers.h>
#include <Math.d/Vector.h>
#include <Math.d/CuCSparse.h>
#include <Timers.d/GetTime.h>
#include <Control.d/ControlInterface.h>
#include <Threads.d/PHelper.h>
#include <Paral.d/SubDOp.h>
#ifdef DISTRIBUTED
#include <Utils.d/DistHelper.h>
#include <Dist.d/DistDom.h>
#endif
#include <Driver.d/DecDomain.h>
#include <Hetero.d/DistFlExchange.h>
#include <Utils.d/SolverInfo.h>
#include <Feti.d/Feti.h>
#include <Utils.d/ModeData.h>

#define EXPLICIT_UPDATE

extern ModeData modeData;

MultiDomainOp::MultiDomainOp(void (MultiDomainOp::*_f)(int),  SubDomain **_sd,
                             DistrVector *_v1, DistrVector*_v2, double c,
                             double *_userDefDisps, SubDOp* _Kuc)
{
 v1 = _v1;
 v2 = _v2;
 f  = _f;
 c1 = c;
 sd = _sd;
 userDefDisps = _userDefDisps;
 Kuc = _Kuc;
}

MultiDomainOp::MultiDomainOp(void (MultiDomainOp::*_f)(int),  SubDomain **_sd,
                             DistrVector *_v1, SubDOp* _Kuc)
{
 v1 = _v1;
 f  = _f;
 sd = _sd;
 Kuc = _Kuc;
}

MultiDomainOp::MultiDomainOp(void (MultiDomainOp::*_f)(int),  SubDomain **_sd)
{
  f = _f;
 sd = _sd;
}

MultiDomainOp::MultiDomainOp(void (MultiDomainOp::*_f)(int),  SubDomain **_sd,
                             DistrVector *_v1, DistrVector*_v2, DistrVector *_v3)
{
 v1 = _v1;
 v2 = _v2;
 v3 = _v3;
 f  = _f;
 sd = _sd;
}

MultiDomainOp::MultiDomainOp(void (MultiDomainOp::*_f)(int),  SubDomain **_sd,
                             DistrVector *_v1, DistrVector*_v2, DistrVector *_v3,  DistrVector *_v4)
{
 v1 = _v1;
 v2 = _v2;
 v3 = _v3;
 v4 = _v4;
 f  = _f;
 sd = _sd;
}

MultiDomainOp::MultiDomainOp(void (MultiDomainOp::*_f)(int),  SubDomain **_sd,
                             DistrVector *_v1, DistrVector *_v2)
{
 v1 = _v1;
 v2 = _v2;
 f  = _f;
 sd = _sd;
}


void
MultiDomainOp::runFor(int isub)
{
  (this->*f)(isub);
}

void
MultiDomainOp::computeExtForce(int isub)
{
  // Get the pointer to the part of the vector f corresponding to subdomain isub
  StackVector localf(v1->subData(isub),v1->subLen(isub));

  // Get the pointer to the part of the vector cnst_f corresponding to subdomain isub
  StackVector localg(v2->subData(isub),v2->subLen(isub));

  int *userDataMap = sd[isub]->getUserDispDataMap();
  sd[isub]->computeExtForce4(localf, localg, c1, (*Kuc)[isub]/*, userDefDisps, userDataMap*/); //XXXX broken
}

void
MultiDomainOp::getConstForce(int isub)
{
  // Get the pointer to the part of the vector f correspoding to subdomain sNum
  StackVector f(v1->subData(isub),v1->subLen(isub));
  sd[isub]->computeConstantForce(f, (*Kuc)[isub]);
}

void
MultiDomainOp::makeAllDOFs(int isub)
{
  sd[isub]->makeAllDOFs();
}

void
MultiDomainOp::getInitState(int isub)
{
 StackVector disp( v1->subData(isub),v1->subLen(isub));
 StackVector veloc(v2->subData(isub),v1->subLen(isub));
 StackVector accel(v3->subData(isub),v1->subLen(isub));
 StackVector   v_p(v4->subData(isub),v1->subLen(isub));

 sd[isub]->initDispVeloc(disp, veloc, accel, v_p);
}

void
MultiDomDynPostProcessor::setPostProcessor(DistFlExchanger *exchanger) 
{
  distFlExchanger = exchanger;
}

void
MultiDomDynPostProcessor::setUserDefs(double **disps, double **vels) 
{
  usrDefDisps = disps;
  usrDefVels = vels;
}

void
MultiDomDynPostProcessor::setNodalTemps(DistrVector* _nodalTemps)
{
  nodalTemps = _nodalTemps;
}

void
MultiDomDynPostProcessor::dynamOutput(int tIndex, MDDynamMat &dynOps, DistrVector &distForce, 
                                      DistrVector *distAeroF, SysState<DistrVector>& distState)
{
  if(!times) times = new StaticTimers;
  startTimerMemory(times->output, times->memoryOutput);

  // PJSA 4-15-08: update bcx for time dependent prescribed displacements and velocities (previously done in computeExternalForce2)
  ControlLawInfo *claw = geoSource->getControlLaw();
  ControlInterface *userSupFunc = domain->getUserSuppliedFunction();
  if(claw && claw->numUserDisp) {
    double t = double(tIndex)*domain->solInfo().getTimeStep();
    double *userDefineDisp = new double[claw->numUserDisp];
    double *userDefineVel  = new double[claw->numUserDisp];
    //cerr << "getting usdd at time " << t << " for dynamOutput\n";
    userSupFunc->usd_disp(t,userDefineDisp,userDefineVel);
    paralApply(decDomain->getNumSub(), decDomain->getAllSubDomains(), &GenSubDomain<double>::setUserDefBC, userDefineDisp, userDefineVel);
    delete [] userDefineDisp; delete [] userDefineVel;
  }

  // Send displacements to fluid code (except C0)
  SolverInfo& sinfo = domain->solInfo();
  if(sinfo.aeroFlag >= 0 && !sinfo.lastIt && tIndex != sinfo.initialTimeIndex && sinfo.aeroFlag != 20) {
    // Send u + IDISP6 to fluid code.
    // IDISP6 is used to compute pre-stress effects.
    DistrVector d_n_aero(distState.getDisp());

    if(domain->solInfo().gepsFlg == 1) {

      for(int i = 0; i < decDomain->getNumSub(); ++i) {
        SubDomain *sd = decDomain->getSubDomain(i);
        BCond* iDis6 = sd->getInitDisp6();
        for(int j = 0; j < sd->numInitDisp6(); ++j) {
          int dof = sd->getCDSA()->locate(iDis6[j].nnum, 1 << iDis6[j].dofnum);
          if(dof >= 0)
            d_n_aero.subData(i)[dof] += iDis6[j].val;
        }
      }
    }

    SysState<DistrVector> state(d_n_aero, distState.getVeloc(), distState.getAccel(), distState.getPrevVeloc());

    distFlExchanger->sendDisplacements(state, usrDefDisps, usrDefVels);
    if(verboseFlag) filePrint(stderr, " ... [E] Sent displacements ...\n");
  }

  if(sinfo.aeroheatFlag >= 0 && tIndex != 0) {
    SysState<DistrVector> tempState(distState.getDisp(), distState.getVeloc(), distState.getPrevVeloc());

    distFlExchanger->sendTemperature(tempState);
    if(verboseFlag) filePrint(stderr, " ... [T] Sent temperatures (%e) ...\n", distState.getDisp().sqNorm());
  }

  if(sinfo.thermohFlag >= 0 && tIndex != 0) {

    for(int i = 0; i < decDomain->getNumSub(); ++i) {
      SubDomain *sd = decDomain->getSubDomain(i);
      for(int j = 0; j < sd->numNodes(); ++j) {
        int tloc  = sd->getCDSA()->locate(j, DofSet::Temp);
        int tloc1 = sd->getDSA()->locate(j, DofSet::Temp);
        double temp  = (tloc >= 0) ? distState.getDisp().subData(i)[tloc] : sd->getBcx()[tloc1];
        if(tloc1 < 0) temp = 0.0;
        nodalTemps->subData(i)[j] = temp;
      }
    }

    distFlExchanger->sendStrucTemp(*nodalTemps);
    if(verboseFlag) filePrint(stderr," ... [T] Sent temperatures ...\n");
  }

  decDomain->postProcessing(distState.getDisp(), distForce, 0.0, distAeroF, tIndex, &dynOps, &distState); 
  stopTimerMemory(times->output, times->memoryOutput);

//  SolverInfo& sinfo = domain->solInfo();
//  if(sinfo.aeroFlag >= 0)
//    if(tIndex != sinfo.initialTimeIndex)
//      distFlExchanger->sendDisplacements(distState, usrDefDisps, usrDefVels);
}

MultiDomainDynam::~MultiDomainDynam() 
{ 
  int nsub = decDomain->getNumSub(); 
  delete decDomain; 
  delete times; 
  if(geomState) delete geomState; 
  if(kelArray) { for(int i=0; i<nsub; ++i) delete [] kelArray[i]; delete [] kelArray; }
  if(allCorot) { for(int i=0; i<nsub; ++i) delete [] allCorot[i]; delete [] allCorot; } 
  if(dprev) delete dprev;
}

MDDynamMat *
MultiDomainDynam::buildOps(double coeM, double coeC, double coeK)
{
  // Have each subdomain create their operators, then put the
  // dynamic matrices in the Feti Solver
  dynMat = new MDDynamMat;

  times->getFetiSolverTime -= getTime(); // PJSA 5-25-05
  decDomain->buildOps(*dynMat, coeM, coeC, coeK, (Rbm **) 0, kelArray);

  if(domain->tdenforceFlag()) { 
    domain->MakeNodalMass(dynMat->M, decDomain->getAllSubDomains());
  }

  times->getFetiSolverTime += getTime();
  return dynMat;
}

MultiDomainDynam::MultiDomainDynam(Domain *d)
{
  domain = d;

#ifdef DISTRIBUTED
  decDomain = new GenDistrDomain<double>(domain);
#else
  decDomain = new GenDecDomain<double>(domain);
#endif
  times  = new StaticTimers;

  claw = 0;
  userSupFunc = 0;
  kelArray = 0;
  allCorot = 0;
  geomState = 0;
  dprev = 0;
}
                                                                                                 
const DistrInfo &
MultiDomainDynam::solVecInfo()
{
  return decDomain->solVecInfo();
}
                                                                                                 
DistrInfo &
MultiDomainDynam::bcInfo()
{
  // prescribed boundary condition distributed vector information
  return *decDomain->pbcVectorInfo();
}

void
MultiDomainDynam::processLastOutput()
{
  OutputInfo *oinfo = geoSource->getOutputInfo();
  for (int iOut = 0; iOut < geoSource->getNumOutInfo(); iOut++)
    oinfo[iOut].interval = 1;
}
                                                                                                 
void
MultiDomainDynam::preProcess()
{
  times->preProcess -= getTime();
                                                                                                 
  // Makes local renumbering, connectivities and dofsets
  decDomain->preProcess();

  // Make all element's dofs
  MultiDomainOp mdop(&MultiDomainOp::makeAllDOFs, decDomain->getAllSubDomains());
#ifdef DISTRIBUTED
  execParal(decDomain->getNumSub(), &mdop, &MultiDomainOp::runFor);
#else
  threadManager->execParal(decDomain->getNumSub(), &mdop);
#endif
  times->preProcess += getTime();

  // Check for user supplied routines (control, force or displacement)
  claw = geoSource->getControlLaw();
  userSupFunc = geoSource->getUserSuppliedFunction();

  // Make the geomState (used for prestress, explicit geometric nonlinear and contact)
  if((domain->solInfo().gepsFlg == 1 && domain->numInitDisp6() > 0) || domain->solInfo().isNonLin() || domain->tdenforceFlag()) {
    times->timeGeom -= getTime();
    geomState = new DistrGeomState(decDomain);
    times->timeGeom += getTime();
  }

  // Update geomState with prescribed dirichlet boundary conditions (explicit geometric nonlinear and contact)
  if(domain->solInfo().isNonLin() || domain->tdenforceFlag())
    execParal(decDomain->getNumSub(), this, &MultiDomainDynam::initSubPrescribedDisplacement);

  // Make corotators and kelArray (used for prestress and explicit geometric nonlinear)
  if((domain->solInfo().gepsFlg == 1 && domain->numInitDisp6() > 0) || domain->solInfo().isNonLin()) {
    times->corotatorTime -= getTime();
    allCorot = new Corotator**[decDomain->getNumSub()];
    execParal(decDomain->getNumSub(), this, &MultiDomainDynam::makeSubCorotators);
    times->corotatorTime += getTime();

    kelArray = new FullSquareMatrix*[decDomain->getNumSub()];
    execParal(decDomain->getNumSub(), this, &MultiDomainDynam::makeSubElementArrays);
  }

  // Initialization for contact
  if(domain->tdenforceFlag())
    domain->InitializeDynamicContactSearch(decDomain->getNumSub(), decDomain->getAllSubDomains());
/* currently done in main.C for linear case
  else
    domain->InitializeStaticContactSearch(MortarHandler::CTC, decDomain->getNumSub(), decDomain->getAllSubDomains());
*/
}

void
MultiDomainDynam::makeSubCorotators(int isub)
{
  SubDomain *sd  = decDomain->getSubDomain(isub);
  int numele     = sd->numElements();
  allCorot[isub] = new Corotator*[numele];
  sd->createCorotators(allCorot[isub]);
}

void
MultiDomainDynam::makeSubElementArrays(int isub)
{
  SubDomain *sd = decDomain->getSubDomain(isub);

  // allocate the element stiffness array
  sd->createKelArray(kelArray[isub]);
 
  // update geomState with IDISP6 if GEPS is requested (geometric prestress / linear only)
  if((sd->numInitDisp6() > 0) && (domain->solInfo().gepsFlg == 1)) // GEPS
   (*geomState)[isub]->updatePrescribedDisplacement(sd->getInitDisp6(), sd->numInitDisp6());

  // build the element stiffness matrices.
  Vector elementInternalForce(sd->maxNumDOF(), 0.0);
  Vector residual(sd->numUncon(), 0.0);
  sd->getStiffAndForce(*(*geomState)[isub], elementInternalForce, allCorot[isub], kelArray[isub], residual);
}

void
MultiDomainDynam::initSubPrescribedDisplacement(int isub)
{
  SubDomain *sd = decDomain->getSubDomain(isub);
  if(sd->nDirichlet() > 0) 
    (*geomState)[isub]->updatePrescribedDisplacement(sd->getDBC(), sd->nDirichlet());
}

void
MultiDomainDynam::getTimes(double &dt, double &tmax)
{
  dt   = domain->solInfo().getTimeStep();
  tmax = domain->solInfo().tmax;
}

void
MultiDomainDynam::getInitialTime(int &initTimeIndex, double &initTime)
{
  initTimeIndex = domain->solInfo().initialTimeIndex;
  initTime      = domain->solInfo().initialTime;
}

double
MultiDomainDynam::getInitialForceNorm()
{
  return domain->solInfo().initExtForceNorm;
}

int
MultiDomainDynam::getTimeIntegration()
{
  return domain->solInfo().timeIntegration;
}

#include <algorithm>
int
MultiDomainDynam::getFilterFlag()
{
  return std::max(domain->solInfo().hzemFilterFlag, domain->solInfo().filterFlags);
}

int*
MultiDomainDynam::boundary()
{
  filePrint(stderr, "Paral.d/MDDynam.C: boundary not implemented here\n");
  return 0;
}

double*
MultiDomainDynam::boundaryValue()
{
  filePrint(stderr, "Paral.d/MDDynam.C: boundaryValue not implemented here\n");
  return 0;
}

Domain*
MultiDomainDynam::getDomain()
{
  return domain;
}
                                                                                                 
void
MultiDomainDynam::getNewMarkParameters(double &beta, double &gamma,
                                       double &alphaf, double &alpham)
{
  beta  = domain->solInfo().newmarkBeta;
  gamma = domain->solInfo().newmarkGamma;
  alphaf = domain->solInfo().newmarkAlphaF;
  alpham = domain->solInfo().newmarkAlphaM;
}
                                                                                                 
void
MultiDomainDynam::getQuasiStaticParameters(double &maxVel, double &delta)
{
  maxVel  = domain->solInfo().qsMaxvel;
  delta  = domain->solInfo().delta;
}
                                                                                                 
void
MultiDomainDynam::getSteadyStateParam(int &steadyFlag, int &steadyMin,
                                      int &steadyMax, double &steadyTol)
{
  steadyFlag  = domain->solInfo().steadyFlag;
  steadyMin   = domain->solInfo().steadyMin;
  steadyMax   = domain->solInfo().steadyMax;
  steadyTol   = domain->solInfo().steadyTol;
}

void
MultiDomainDynam::getContactForce(DistrVector &d, DistrVector &ctc_f)
{
  // DEBUG CONTACT
  times->tdenforceTime -= getTime();
  ctc_f.zero();
  if(domain->tdenforceFlag()) {
    times->updateSurfsTime -= getTime();
    domain->UpdateSurfaces(geomState, 1, decDomain->getAllSubDomains()); // update to current configuration
    times->updateSurfsTime += getTime();
#ifdef EXPLICIT_UPDATE
    execParal1R(decDomain->getNumSub(), this, &MultiDomainDynam::subExplicitUpdate, d);
#else
    DistrVector dinc(decDomain->solVecInfo());
    dinc.linC(1.0, d, -1.0, *dprev);
    geomState->update(dinc);
    (*dprev) = d;
#endif
    times->updateSurfsTime -= getTime();
    domain->UpdateSurfaces(geomState, 2, decDomain->getAllSubDomains()); // update to predicted configuration
    times->updateSurfsTime += getTime();

    times->contactSearchTime -= getTime();
    domain->PerformDynamicContactSearch(domain->solInfo().getTimeStep());
    times->contactSearchTime += getTime();

    times->contactForcesTime -= getTime();
    domain->AddContactForces(domain->solInfo().getTimeStep(), ctc_f);
    times->contactForcesTime += getTime();
  }
  times->tdenforceTime += getTime();
}

void
MultiDomainDynam::computeExtForce2(SysState<DistrVector> &distState,
                                   DistrVector &f, DistrVector &cnst_f, int tIndex,
                                   double t, DistrVector *aero_f,
                                   double gamma, double alphaf)
{
  times->formRhs -= getTime();

  // compute USDD prescribed displacements
  // I think we should call the control functions on the rank zero mpi process then send to all others
  double *userDefineDisp = 0;
  if(claw && userSupFunc) {
    if(claw->numUserDisp) {
      userDefineDisp = new double[claw->numUserDisp];
      double *userDefineVel  = new double[claw->numUserDisp];
      userSupFunc->usd_disp(t, userDefineDisp, userDefineVel);
      paralApply(decDomain->getNumSub(), decDomain->getAllSubDomains(), &GenSubDomain<double>::setUserDefBC, userDefineDisp, userDefineVel); // update bcx, vcx
      delete [] userDefineVel;
    }
  }

  // update geomState for nonlinear problems. note computeExtForce2 must be called before getInternalForce
  if(domain->solInfo().isNonLin() || domain->tdenforceFlag()) {
#ifdef EXPLICIT_UPDATE
    execParal1R(decDomain->getNumSub(), this, &MultiDomainDynam::subExplicitUpdate, distState.getDisp());
#else
    if(!dprev) { dprev = new DistrVector(decDomain->solVecInfo()); dprev->zero(); }
    DistrVector dinc(decDomain->solVecInfo());
    dinc.linC(1.0, distState.getDisp(), -1.0, *dprev); // incremental displacement: dinc = d - dprev
    geomState->update(dinc);
    *dprev = distState.getDisp();
#endif
    execParal1R(decDomain->getNumSub(), this, &MultiDomainDynam::subUpdateGeomStateUSDD, userDefineDisp);
    geomState->setVelocity(distState.getDisp(), distState.getVeloc(), distState.getAccel());
  }

  // update nodal temperatures for thermoe problem
  if(domain->solInfo().thermoeFlag >= 0 && tIndex >= 0) {
    distFlExchanger->getStrucTemp(nodalTemps->data());
    if(verboseFlag) filePrint(stderr," ... [E] Received temperatures ...\n");
  }


  // add f(t) to cnst_f
  MultiDomainOp mdop(&MultiDomainOp::computeExtForce,
                     decDomain->getAllSubDomains(), &f, &cnst_f, t, userDefineDisp, dynMat->Kuc);
  threadManager->execParal(decDomain->getNumSub(), &mdop);
  if(userDefineDisp) delete [] userDefineDisp;

  // add USDF forces
  if(claw && userSupFunc) {
    if(claw->numUserForce) {
      double *userDefineForce = new double[claw->numUserForce];
      userSupFunc->usd_forc(t, userDefineForce);
      decDomain->addUserForce(f, userDefineForce);
      delete [] userDefineForce;
    }
  }

  // add ACTUATOR forces XXXX but is distState at time t? need to check disp/vel/acc for explicit/implicit/quasistatics
  if(claw && userSupFunc) {
    if(claw->numActuator) {
      double *ctrdisp = new double[claw->numSensor];
      double *ctrvel  = new double[claw->numSensor];
      double *ctracc  = new double[claw->numSensor];
      double *ctrfrc  = new double[claw->numActuator];
#ifdef DISTRIBUTED
      for(int i=0; i<claw->numSensor; ++i) ctrdisp[i] = ctrvel[i] = ctracc[i] = std::numeric_limits<double>::min();
#endif
      DistrVector &disp = distState.getDisp();
      DistrVector &vel = distState.getVeloc();
      DistrVector &acc = distState.getAccel();
      decDomain->extractControlData(disp, vel, acc, ctrdisp, ctrvel, ctracc);
#ifdef DISTRIBUTED
      structCom->globalMax(claw->numSensor, ctrdisp);
      structCom->globalMax(claw->numSensor, ctrvel);
      structCom->globalMax(claw->numSensor, ctracc);
#endif
      userSupFunc->ctrl(ctrdisp, ctrvel, ctracc, ctrfrc, t);
      decDomain->addCtrl(f, ctrfrc);
      delete [] ctrdisp; delete [] ctrvel; delete [] ctracc; delete [] ctrfrc;
    }
  }

  // add aeroelastic forces from fluid dynamics code
  SolverInfo& sinfo = domain->solInfo();
  if(sinfo.aeroFlag >= 0 && tIndex >= 0) {

    aeroForce->zero();
    int iscollocated;
    double tFluid = distFlExchanger->getFluidLoad(*aeroForce, tIndex, t,
                                                  alphaf, iscollocated);
    if(verboseFlag) filePrint(stderr," ... [E] Received fluid forces ...\n");

    if(sinfo.aeroFlag == 20) {
      if(prevIndex >= 0)
        aero_f->linC(0.5,*aeroForce,0.5,*prevFrc);
      else
        *aero_f = *aeroForce;
    }
    else {
      if (iscollocated == 0) {
        if(prevIndex >= 0) {
          *aeroForce *= (1/gamma);
          aeroForce->linAdd(((gamma-1.0)/gamma),*prevFrc);
        }
      }

      double alpha = 1.0-alphaf;
      if(prevIndex < 0) alpha = 1.0;

      aero_f->linC(alpha, *aeroForce, (1.0-alpha), *prevFrc);
    }
    f += *aero_f;

    *prevFrc = *aeroForce;
    prevTime = tFluid;
    prevIndex = tIndex;
  }

  // add aerothermal fluxes from fluid dynamics code
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

  // KHP: apply projector here
  if(domain->solInfo().filterFlags || domain->solInfo().hzemFilterFlag)
    trProject(f);

  if(tIndex == 1)
    domain->solInfo().initExtForceNorm = f.norm();

  times->formRhs += getTime();
}

void
MultiDomainDynam::getConstForce(DistrVector& v)
{
  times->formRhs -= getTime();
  MultiDomainOp mdop(&MultiDomainOp::getConstForce, decDomain->getAllSubDomains(), &v, dynMat->Kuc);
  threadManager->execParal(decDomain->getNumSub(), &mdop);
  times->formRhs += getTime();
}

void
MultiDomainDynam::getInitState(SysState<DistrVector>& state)
{
  // initialize state with IDISP/IDISP6/IVEL/IACC or RESTART (XXXX initial accelerations are currently not supported)
  MultiDomainOp mdop(&MultiDomainOp::getInitState, decDomain->getAllSubDomains(),
                     &state.getDisp(), &state.getVeloc(), &state.getAccel(),
                     &state.getPrevVeloc());
  threadManager->execParal(decDomain->getNumSub(), &mdop);

  // if we have a user supplied function, give it the initial state at the sensors
  // .. first update bcx, vcx in case any of the sensors have prescribed displacements
  if(claw && userSupFunc) {
    if(claw->numUserDisp) {
      double *userDefineDisp = new double[claw->numUserDisp];
      double *userDefineVel  = new double[claw->numUserDisp];
      userSupFunc->usd_disp(domain->solInfo().initialTime, userDefineDisp, userDefineVel);
      paralApply(decDomain->getNumSub(), decDomain->getAllSubDomains(), &GenSubDomain<double>::setUserDefBC, userDefineDisp, userDefineVel);
      delete [] userDefineDisp; delete [] userDefineVel;
    }
    if(claw->numSensor) {
      double *ctrdisp = new double[claw->numSensor];
      double *ctrvel  = new double[claw->numSensor];
      double *ctracc  = new double[claw->numSensor];
#ifdef DISTRIBUTED
      for(int i=0; i<claw->numSensor; ++i) ctrdisp[i] = ctrvel[i] = ctracc[i] = std::numeric_limits<double>::min();
#endif
      DistrVector &disp = state.getDisp();
      DistrVector &vel  = state.getVeloc();
      DistrVector &acc  = state.getAccel();
      decDomain->extractControlData(disp, vel, acc, ctrdisp, ctrvel, ctracc);
#ifdef DISTRIBUTED
      structCom->globalMax(claw->numSensor, ctrdisp);
      structCom->globalMax(claw->numSensor, ctrvel);
      structCom->globalMax(claw->numSensor, ctracc);
#endif
      userSupFunc->init(ctrdisp, ctrvel, ctracc);
      delete [] ctrdisp; delete [] ctrvel; delete [] ctracc;
    }
  }
}

MultiDomDynPostProcessor *
MultiDomainDynam::getPostProcessor()
{
 if(domain->solInfo().aeroFlag >= 0) {
   mddPostPro = new MultiDomDynPostProcessor(decDomain, distFlExchanger, times); 
   return mddPostPro;
 }
 else {
   mddPostPro = new MultiDomDynPostProcessor(decDomain, times);
 }
 return mddPostPro;
}

void
MultiDomainDynam::printTimers(MDDynamMat *dynOps, double timeLoop)
{
  times->numSubdomain = decDomain->getNumSub();
  filePrint(stderr," ... Print Timers                   ... \n");

  if(domain->solInfo().type == 2 && domain->solInfo().fetiInfo.version == 3) {
    times->printFetiDPtimers(domain->getTimers(),
                             dynOps->dynMat->getSolutionTime(),
                             domain->solInfo() ,
                             dynOps->dynMat->getTimers(),
                             geoSource->getCheckFileInfo()[0],
                             domain);
  }
  else {
    times->printStaticTimers(domain->getTimers(),
                             dynOps->dynMat->getSolutionTime(),
                             domain->solInfo() ,
                             dynOps->dynMat->getTimers(),
                             geoSource->getCheckFileInfo()[0],
                             domain);
 }

/*
   switch(domain->solInfo().fetiInfo.version) {
     default:
     case FetiInfo::feti1:
     case FetiInfo::feti2:
       times->printStaticTimers(domain->getTimers(),
                                dynOps->dynMat->getSolutionTime(),
                                domain->solInfo() ,
                                dynOps->dynMat->getTimers(),
                                geoSource->getCheckFileInfo()[0],
                                domain);
       break;
                                                                                                 
     case FetiInfo::fetidp:
       times->printFetiDPtimers(domain->getTimers(),
                                dynOps->dynMat->getSolutionTime(),
                                domain->solInfo() ,
                                dynOps->dynMat->getTimers(),
                                geoSource->getCheckFileInfo()[0],
                                domain);
       break;
   }
*/
}

void
MultiDomainDynam::trProject(DistrVector &)
{
  filePrint(stderr, "Paral.d/MDDynam.C: trProject not implemented here\n");
}

void
MultiDomainDynam::project(DistrVector &)
{
  filePrint(stderr, "Paral.d/MDDynam.C: project not implemented here\n");
}

void
MultiDomainDynam::getRayleighCoef(double& alpha)
{
  alpha = domain->solInfo().alphaDamp;
}

void
MultiDomainDynam::addPrescContrib(SubDOp*, SubDOp*, DistrVector&,
                                  DistrVector&, DistrVector&, DistrVector&,
                                  double time)
{
  filePrint(stderr, "Paral.d/MDDynam.C: addPrescContrib not implemented here\n");
}

SubDOp*
MultiDomainDynam::getpK(MDDynamMat* dynOps)
{
  return dynOps->K;
}

SubDOp*
MultiDomainDynam::getpM(MDDynamMat* dynOps)
{
  return dynOps->M;
}

SubDOp*
MultiDomainDynam::getpC(MDDynamMat* dynOps)
{
  return dynOps->C;
}

void
MultiDomainDynam::computeStabilityTimeStep(double dt, MDDynamMat &dMat)
{
  double dt_c; //time step computed from the matrix K
  dt_c = decDomain->computeStabilityTimeStep(dMat);

  filePrint(stderr," **************************************\n");
  if (domain->solInfo().modifiedWaveEquation) {
    dt_c = 1.73205*dt_c;
    filePrint(stderr," CONDITIONALLY STABLE MODIFIED WAVE EQUATION \n");
  }
  else
    filePrint(stderr," CONDITIONALLY STABLE NEWMARK ALGORITHM \n");
  filePrint(stderr," --------------------------------------\n");
  filePrint(stderr," Specified time step      = %10.4e\n",dt);
  filePrint(stderr," Stability max. time step = %10.4e\n",dt_c);
  filePrint(stderr," **************************************\n");
  if( dt_c < dt ) {
    dt = dt_c;
    filePrint(stderr," Stability max. time step is selected\n");
  } else
    filePrint(stderr," Specified time step is selected\n");
  filePrint(stderr," **************************************\n");

  domain->solInfo().setTimeStep(dt);
}

int
MultiDomainDynam::getModeDecompFlag()
{
  return domain->solInfo().modeDecompFlag;
}

void
MultiDomainDynam::modeDecompPreProcess(SparseMatrix* M)
{
  filePrint(stderr, "Paral.d/MDDynam.C: modeDecompPreProcess not implemented here\n");
}

void
MultiDomainDynam::modeDecomp(double t, int tIndex, DistrVector& d_n)
{
  filePrint(stderr, "Paral.d/MDDynam.C: modeDecomp not implemented here\n");
}

void 
MultiDomainDynam::getInternalForce(DistrVector &d, DistrVector &f, double t)
{
  if(domain->solInfo().isNonLin())  // PJSA 3-31-08
    execParal2R(decDomain->getNumSub(), this, &MultiDomainDynam::subGetInternalForce, f, t);
  else {
    f.zero();
    execParal2R(decDomain->getNumSub(), this, &MultiDomainDynam::subGetKtimesU, d, f);
  }

  if(domain->solInfo().filterFlags || domain->solInfo().hzemFilterFlag)
    trProject(f);
}

void
MultiDomainDynam::subExplicitUpdate(int isub, DistrVector &d)
{
  SubDomain *sd = decDomain->getSubDomain(isub);
  StackVector subd(d.subData(isub), d.subLen(isub));
  (*geomState)[isub]->explicitUpdate(sd->getNodes(), subd);
}


void
MultiDomainDynam::subUpdateGeomStateUSDD(int isub, double *userDefineDisp)
{
  SubDomain *sd = decDomain->getSubDomain(isub);
  ControlLawInfo *subClaw = sd->getClaw();
  if(subClaw) {
    if(subClaw->numUserDisp) {
      double *subUserDefineDisp = new double[subClaw->numUserDisp];
      for(int i=0; i<subClaw->numUserDisp; ++i) 
        subUserDefineDisp[i] = userDefineDisp[sd->getUserDispDataMap()[i]];
      (*geomState)[isub]->updatePrescribedDisplacement(subUserDefineDisp, subClaw, sd->getNodes());
      delete [] subUserDefineDisp;
    }
  } 
}

void
MultiDomainDynam::subGetInternalForce(int isub, DistrVector &f, double t)
{
  SubDomain *sd = decDomain->getSubDomain(isub);
  Vector residual(f.subLen(isub), 0.0);
  Vector eIF(sd->maxNumDOF()); // eIF = element internal force for one element (a working array)
  sd->getStiffAndForce(*(*geomState)[isub], eIF, allCorot[isub], kelArray[isub], residual, 1.0, t); // residual -= internal force
  StackVector subf(f.subData(isub), f.subLen(isub));
  subf.linC(residual,-1.0); // f = -residual
}

void
MultiDomainDynam::subGetKtimesU(int isub, DistrVector &d, DistrVector &f)
{
  SubDomain *sd = decDomain->getSubDomain(isub);
  StackVector subf(f.subData(isub), f.subLen(isub));
  StackVector subd(d.subData(isub), d.subLen(isub));
  sd->getKtimesU(subd, (double *) 0, subf, 1.0, (kelArray) ? kelArray[isub] : (FullSquareMatrix *) 0);
}

void
MultiDomainDynam::computeTimeInfo()
{
  // Time integration information

  // Get total time and time step size and store them
  double totalTime = domain->solInfo().tmax;
  double dt        = domain->solInfo().getTimeStep();

  // Compute maximum number of steps
  int maxStep = (int) ( (totalTime+0.49*dt)/dt );

  // Compute time remainder
  double remainder = totalTime - maxStep*dt;
  if(std::abs(remainder)>0.01*dt){
    domain->solInfo().tmax = maxStep*dt;
    fprintf(stderr, " Warning: Total time is being changed to : %e\n", domain->solInfo().tmax);
  }
}

int
MultiDomainDynam::aeroPreProcess(DistrVector &disp, DistrVector &vel,
                                 DistrVector &accel, DistrVector &lastVel)  
{
  // get solver info
  SolverInfo& sinfo = domain->solInfo();

  // Initialize previous force data
  // Reexamine for the case of restart
  prevFrc = new DistrVector(solVecInfo());
  prevFrcBackup = new DistrVector(solVecInfo());
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

  mddPostPro->setPostProcessor(distFlExchanger);
  mddPostPro->setUserDefs(usrDefDisps, usrDefVels);

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

int 
MultiDomainDynam::cmdCom(int cmdFlag) 
{
  return distFlExchanger->cmdCom(cmdFlag);
}

int
MultiDomainDynam::getAeroAlg()
{
  return domain->solInfo().aeroFlag;
}

void
MultiDomainDynam::aeroSend(double time, DistrVector& d_n, DistrVector& v_n, DistrVector& a_n, DistrVector& v_p)
{
  // Send u + IDISP6 to fluid code.
  // IDISP6 is used to compute pre-stress effects.
  DistrVector d_n_aero(d_n);

  if(domain->solInfo().gepsFlg == 1) {

    for(int i = 0; i < decDomain->getNumSub(); ++i) {
      SubDomain *sd = decDomain->getSubDomain(i);
      BCond* iDis6 = sd->getInitDisp6();
      for(int j = 0; j < sd->numInitDisp6(); ++j) {
        int dof = sd->getCDSA()->locate(iDis6[j].nnum, 1 << iDis6[j].dofnum);
        if(dof >= 0)
          d_n_aero.subData(i)[dof] += iDis6[j].val;
      }
    }
  }

  SysState<DistrVector> state(d_n_aero, v_n, a_n, v_p);

  distFlExchanger->sendDisplacements(state, usrDefDisps, usrDefVels);
  if(verboseFlag) filePrint(stderr, " ... [E] Sent displacements ...\n");
}

void
MultiDomainDynam::a5TimeLoopCheck(int& parity, double& t, double dt)
{
  if(domain->solInfo().aeroFlag == 5) {
     if(!parity) t -= dt;
     parity = ( parity ? 0 : 1 );
  }
}

void
MultiDomainDynam::a5StatusRevise(int parity, SysState<DistrVector>& curState, SysState<DistrVector>& bkState)
{
  if(domain->solInfo().aeroFlag == 5) {
    if(parity) { // restore
      *prevFrc = *prevFrcBackup;
      prevIndex = prevIndexBackup;
      prevTime = prevTimeBackup;
      curState.getDisp()      = bkState.getDisp();
      curState.getVeloc()     = bkState.getVeloc();
      curState.getAccel()     = bkState.getAccel();
      curState.getPrevVeloc() = bkState.getPrevVeloc();
    }
    else { // backup
      *prevFrcBackup = *prevFrc;
      prevIndexBackup = prevIndex;
      prevTimeBackup = prevTime;
      bkState.getDisp()      = curState.getDisp();
      bkState.getVeloc()     = curState.getVeloc();
      bkState.getAccel()     = curState.getAccel();
      bkState.getPrevVeloc() = curState.getPrevVeloc();
    }
  }
}

void
MultiDomainDynam::thermoePreProcess(DistrVector&, DistrVector&, DistrVector&)
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
      mddPostPro->setPostProcessor(distFlExchanger);
      mddPostPro->setUserDefs(usrDefDisps, usrDefVels);
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
MultiDomainDynam::thermohPreProcess(DistrVector& d, DistrVector&, DistrVector&)
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
      mddPostPro->setPostProcessor(distFlExchanger);
      mddPostPro->setUserDefs(usrDefDisps, usrDefVels);
    }

    nodalTemps = new DistrVector(decDomain->ndVecInfo());
    int buffLen = nodalTemps->size();
    mddPostPro->setNodalTemps(nodalTemps);

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

int
MultiDomainDynam::getThermoeFlag()
{
  return domain->solInfo().thermoeFlag;
}

int
MultiDomainDynam::getThermohFlag()
{
  return domain->solInfo().thermohFlag;
}

void
MultiDomainDynam::aeroHeatPreProcess(DistrVector& disp, DistrVector& vel, DistrVector& lastVel)
{
  // get solver info
  SolverInfo& sinfo = domain->solInfo();

  // Initialize previous force data
  // Reexamine for the case of restart
  prevFrc = new DistrVector(solVecInfo());
  prevFrcBackup = new DistrVector(solVecInfo());
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

  mddPostPro->setPostProcessor(distFlExchanger);
  mddPostPro->setUserDefs(usrDefDisps, usrDefVels);

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
MultiDomainDynam::getAeroheatFlag()
{
  return domain->solInfo().aeroheatFlag;
}

