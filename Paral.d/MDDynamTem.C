#include <Paral.d/MDDynam.h>
#include <iostream>
#include <Math.d/Skyline.d/SkyMatrix.h>
#include <Paral.d/MDOp.h>
#include <Timers.d/StaticTimers.h>
#include <Math.d/Vector.h>
#include <Math.d/SparseMatrix.h>
#include <Math.d/CuCSparse.h>
#include <Math.d/NBSparseMatrix.h>
#include <Math.d/DBSparseMatrix.h>
#include <Math.d/SGISparseMatrix.h>
#include <Math.d/BLKSparseMatrix.h>
#include <Math.d/Skyline.d/SGISky.h>
#include <Timers.d/GetTime.h>
#include <Control.d/ControlInterface.h>
#include <Driver.d/DecDomain.h>
#include <Feti.d/Feti.h>
#include <Driver.d/StructProp.h>
#include <Dist.d/DistDom.h>
#include <Hetero.d/DistFlExchange.h>
#include <Comm.d/Communicator.h>
#include <Utils.d/SolverInfo.h>

template <class Scalar>
MultiDomainDynam<Scalar>::~MultiDomainDynam() 
{ 
  int nsub = decDomain->getNumSub(); 
  delete decDomain; 
  delete times; 
  if(geomState) delete geomState; 
  if(kelArray) { for(int i=0; i<nsub; ++i) delete [] kelArray[i]; delete [] kelArray; }
  if(allCorot) { for(int i=0; i<nsub; ++i) delete [] allCorot[i]; delete [] allCorot; } 
  if(dprev) delete dprev;
}

template <class Scalar>
MDDynamMat *
MultiDomainDynam<Scalar>::buildOps(double coeM, double coeC, double coeK)
{
 // Have each subdomain create their operators, then put the
 // dynamic matrices in the Feti Solver
 MDDynamMat *dynMat = new MDDynamMat;

 times->getFetiSolverTime -= getTime(); // PJSA 5-25-05
 decDomain->buildOps(*dynMat, coeM, coeC, coeK, (Rbm **) 0, kelArray);

 if(domain->tdenforceFlag()) { 
   domain->MakeNodalMass(dynMat->M, decDomain->getAllSubDomains());
 }

 times->getFetiSolverTime += getTime();
 return dynMat;
}

template <class Scalar>
MultiDomainDynam<Scalar>::MultiDomainDynam(Domain *d)
{
 domain = d;

 switch(domain->solInfo().type) {
   default:
   case 2: // FETI
     switch(domain->solInfo().fetiInfo.version) {
       default:
       case FetiInfo::feti1:
         filePrint(stderr, " ... FETI-1 is Selected             ...\n");
         break;
       case FetiInfo::feti2:
         filePrint(stderr, " ... FETI-2 is Selected             ...\n");
         break;
       case FetiInfo::fetidp:
         if (!(d->solInfo().fetiInfo.dph_flag))
           filePrint(stderr, " ... FETI-Dual/Primal is Selected   ...\n");
         else
           filePrint(stderr, " ... FETI-DPH is Selected           ...\n");
         break;
     }
     break;
   case 3: // "Block Diagonal" solver
     filePrint(stderr, " ... Diag Parallel is Selected      ...\n");
     break;
  }
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
                                                                                                 
template <class Scalar>
const DistrInfo &
MultiDomainDynam<Scalar>::solVecInfo()
{
 return decDomain->solVecInfo();
}
                                                                                                 
template <class Scalar>
DistrInfo &
MultiDomainDynam<Scalar>::bcInfo()
{
 // prescribed boundary condition distributed vector information
 return *decDomain->pbcVectorInfo();
}

template <class Scalar>
void
MultiDomainDynam<Scalar>::processLastOutput()  {

  OutputInfo *oinfo = geoSource->getOutputInfo();
  for (int iOut = 0; iOut < geoSource->getNumOutInfo(); iOut++)
    oinfo[iOut].interval = 1;
}
                                                                                                 
template <class Scalar>
void
MultiDomainDynam<Scalar>::preProcess()
{
  times->preProcess -= getTime();
                                                                                                 
  // Makes renumbering, connectivities and dofsets
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
    execParal(decDomain->getNumSub(), this, &MultiDomainDynam<Scalar>::initSubPrescribedDisplacement);

  // Make corotators and kelArray (used for prestress and explicit geometric nonlinear)
  if((domain->solInfo().gepsFlg == 1 && domain->numInitDisp6() > 0) || domain->solInfo().isNonLin()) {
    times->corotatorTime -= getTime();
    allCorot = new Corotator**[decDomain->getNumSub()];
    execParal(decDomain->getNumSub(), this, &MultiDomainDynam<Scalar>::makeSubCorotators);
    times->corotatorTime += getTime();

    kelArray = new FullSquareMatrix*[decDomain->getNumSub()];
    execParal(decDomain->getNumSub(), this, &MultiDomainDynam<Scalar>::makeSubElementArrays);
  }

  // Initialization for contact
  if(domain->tdenforceFlag())
    domain->InitializeDynamicContactSearch(decDomain->getNumSub(), decDomain->getAllSubDomains());
}

template <class Scalar>
void
MultiDomainDynam<Scalar>::makeSubCorotators(int isub)
{
  SubDomain *sd  = decDomain->getSubDomain(isub);
  int numele     = sd->numElements();
  allCorot[isub] = new Corotator*[numele];
  sd->createCorotators(allCorot[isub]);
}

template <class Scalar>
void
MultiDomainDynam<Scalar>::makeSubElementArrays(int isub)
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

template <class Scalar>
void
MultiDomainDynam<Scalar>::initSubPrescribedDisplacement(int isub)
{
  SubDomain *sd = decDomain->getSubDomain(isub);
  if(sd->nDirichlet() > 0) 
    (*geomState)[isub]->updatePrescribedDisplacement(sd->getDBC(), sd->nDirichlet());
}

template <class Scalar>
void
MultiDomainDynam<Scalar>::getTimes(double &dt, double &tmax)
{
  dt   = domain->solInfo().dt;
  tmax = domain->solInfo().tmax;
}

template <class Scalar>
void
MultiDomainDynam<Scalar>::getInitialTime(int &initTimeIndex, double &initTime)
{
  initTimeIndex = domain->solInfo().initialTimeIndex;
  initTime      = domain->solInfo().initialTime;
}

template <class Scalar>
double
MultiDomainDynam<Scalar>::getInitialForceNorm()
{
  return domain->solInfo().initExtForceNorm;
}

template <class Scalar>
int
MultiDomainDynam<Scalar>::getTimeIntegration()
{
  return domain->solInfo().timeIntegration;
}
                                                                                                 
template <class Scalar>
int
MultiDomainDynam<Scalar>::getFilterFlag()
{
  return domain->solInfo().filterFlags;
}
                                                                                                 
template <class Scalar>
void
MultiDomainDynam<Scalar>::getNewMarkParameters(double &beta, double &gamma,
                                       double &alphaf, double &alpham)
{
  beta  = domain->solInfo().newmarkBeta;
  gamma = domain->solInfo().newmarkGamma;
  alphaf = domain->solInfo().newmarkAlphaF;
  alpham = domain->solInfo().newmarkAlphaM;
}
                                                                                                 
template <class Scalar>
void
MultiDomainDynam<Scalar>::getQuasiStaticParameters(double &maxVel, double &delta)
{
  maxVel  = domain->solInfo().qsMaxvel;
  delta  = domain->solInfo().delta;
}
                                                                                                 
template <class Scalar>
void
MultiDomainDynam<Scalar>::getSteadyStateParam(int &steadyFlag, int &steadyMin,
                                         int &steadyMax, double &steadyTol)
{
  steadyFlag  = domain->solInfo().steadyFlag;
  steadyMin   = domain->solInfo().steadyMin;
  steadyMax   = domain->solInfo().steadyMax;
  steadyTol   = domain->solInfo().steadyTol;
}

template <class Scalar>
void
MultiDomainDynam<Scalar>::getContactForce(DistrVector &d, DistrVector &ctc_f)
{
  // DEBUG CONTACT
  times->tdenforceTime -= getTime();
  ctc_f.zero();
  if(domain->tdenforceFlag()) {
    times->updateSurfsTime -= getTime();
    domain->UpdateSurfaces(geomState, 1, decDomain->getAllSubDomains()); // update to current configuration
    times->updateSurfsTime += getTime();

    DistrVector dinc(decDomain->solVecInfo());
    dinc.linC(1.0, d, -1.0, *dprev);
    geomState->update(dinc);

    times->updateSurfsTime -= getTime();
    domain->UpdateSurfaces(geomState, 2, decDomain->getAllSubDomains()); // update to predicted configuration
    times->updateSurfsTime += getTime();

    times->contactSearchTime -= getTime();
    domain->PerformDynamicContactSearch(domain->solInfo().dt);
    times->contactSearchTime += getTime();

    times->contactForcesTime -= getTime();
    domain->AddContactForces(domain->solInfo().dt, ctc_f);
    times->contactForcesTime += getTime();

    (*dprev) = d;
  }
  times->tdenforceTime += getTime();
}

template <class Scalar>
void MultiDomainDynam<Scalar>::computeExtForce2(SysState<DistrVector> &distState,
                                                DistrVector &f, DistrVector &cnst_f, int tIndex,
                                                double t, DistrVector *aero_f,
                                                double gamma, double alphaf)
{
  times->formRhs -= getTime();

  // compute USDD prescribed displacements
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
    if(!dprev) { dprev = new DistrVector(decDomain->solVecInfo()); dprev->zero(); }
    DistrVector dinc(decDomain->solVecInfo());
    dinc.linC(1.0, distState.getDisp(), -1.0, *dprev); // incremental displacement: dinc = d - dprev
    geomState->update(dinc);
    execParal1R(decDomain->getNumSub(), this, &MultiDomainDynam<Scalar>::subUpdateGeomStateUSDD, userDefineDisp);
    *dprev = distState.getDisp();
  }

  // add FORCE (including MFTT), HDNB, ROBIN, GRAVITY and PRESSURE forces
  // for linear problems also add non-homogeneous DISP/TEMP and USDD forces
  MultiDomainOp mdop(&MultiDomainOp::computeExtForce,
                     decDomain->getAllSubDomains(), &f, &cnst_f, t, userDefineDisp, geomState);
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

  // this is done in OpMake computeExtForce4 for single domain 
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

  // KHP: apply projector here
  int useProjector = domain->solInfo().filterFlags;
  if(useProjector) trProject(f);

  if(tIndex == 1)
    domain->solInfo().initExtForceNorm = f.norm();

  times->formRhs += getTime();
}

template <class Scalar>
void
MultiDomainDynam<Scalar>::getConstForce(DistrVector& v)
{
  times->formRhs -= getTime();
  MultiDomainOp mdop(&MultiDomainOp::getConstForce, decDomain->getAllSubDomains(), &v);
  threadManager->execParal(decDomain->getNumSub(), &mdop);
  times->formRhs += getTime();
}

template <class Scalar>
void
MultiDomainDynam<Scalar>::getInitState(SysState<DistrVector>&state)
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

template <class Scalar>
MultiDomDynPostProcessor *
MultiDomainDynam<Scalar>::getPostProcessor()
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

template <class Scalar>
void
MultiDomainDynam<Scalar>::printTimers(MDDynamMat *dynOps)
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

template <class Scalar>
void
MultiDomainDynam<Scalar>::addPrescContrib(SubDOp*, SubDOp*, DistrVector&,
                                          DistrVector&, DistrVector&, DistrVector&,
                                          double time)
{
  filePrint(stderr,"Paral.d/MDDynamTem.C: addPrescContrib not implemented here");
}

template <class Scalar>
void
MultiDomainDynam<Scalar>::computeStabilityTimeStep(double dt, MDDynamMat &dMat)
{
  double dt_c;//time step computed from the matrix K
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
  fflush(stderr);
}

template <class Scalar>
void 
MultiDomainDynam<Scalar>::getInternalForce(DistrVector &d, DistrVector &f, double t)
{
  if(domain->solInfo().isNonLin())  // PJSA 3-31-08
    execParal1R(decDomain->getNumSub(), this, &MultiDomainDynam<Scalar>::subGetInternalForce, f);
  else {
    f.zero();
    execParal2R(decDomain->getNumSub(), this, &MultiDomainDynam<Scalar>::subGetKtimesU, d, f);
  }
  int useProjector = domain->solInfo().filterFlags;
  if(useProjector) trProject(f);
}

template <class Scalar>
void
MultiDomainDynam<Scalar>::subUpdateGeomStateUSDD(int isub, double *userDefineDisp)
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

template <class Scalar>
void
MultiDomainDynam<Scalar>::subGetInternalForce(int isub, DistrVector &f)
{
  SubDomain *sd = decDomain->getSubDomain(isub);
  Vector residual(f.subLen(isub), 0.0);
  Vector eIF(sd->maxNumDOF()); // eIF = element internal force for one element (a working array)
  sd->getStiffAndForce(*(*geomState)[isub], eIF, allCorot[isub], kelArray[isub], residual); // residual -= internal force
  StackVector subf(f.subData(isub), f.subLen(isub));
  subf.linC(residual,-1.0); // f = -residual
}

template <class Scalar>
void
MultiDomainDynam<Scalar>::subGetKtimesU(int isub, DistrVector &d, DistrVector &f)
{
  SubDomain *sd = decDomain->getSubDomain(isub);
  StackVector subf(f.subData(isub), f.subLen(isub));
  StackVector subd(d.subData(isub), d.subLen(isub));
  sd->getKtimesU(subd, (double *) 0, subf, 1.0, (kelArray) ? kelArray[isub] : (FullSquareMatrix *) 0);
}

template <class Scalar>
int MultiDomainDynam<Scalar>::aeroPreProcess(DistrVector &disp, DistrVector &vel,
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
  mddPostPro->setPostProcessor(distFlExchanger);
  mddPostPro->setUserDefs(usrDefDisps, usrDefVels);

  // negotiate with the fluid code
  distFlExchanger->negotiate();

  //int restartinc = (solInfo().nRestart >= 0) ? (solInfo().nRestart) : 0;
  // TODO
  int restartinc = 10;

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

  return 0;
}

template <class Scalar>
int MultiDomainDynam<Scalar>::cmdCom(int cmdFlag) 
{
  // return distFlExchanger->cmdCom(cmdFlag);
  fprintf(stderr, "ERROR: Optimization is not available for Feti\n");
  return -1;
}

