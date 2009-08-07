#include <stdlib.h>
#include <stdio.h>
#include <Utils.d/dbg_alloca.h>

#include <Driver.d/Domain.h>
#include <Driver.d/Dynam.h>

#include <Problems.d/DynamDescr.h>

#include <Math.d/FullMatrix.h>
#include <Math.d/SparseMatrix.h>
#include <Math.d/DBSparseMatrix.h>
#include <Math.d/NBSparseMatrix.h>
#include <Math.d/CuCSparse.h>
#include <Math.d/Skyline.d/SkyMatrix.h>

#include <Utils.d/dofset.h>
#include <Solvers.d/Solver.h>
#include <Element.d/State.h>
#include <Timers.d/StaticTimers.h>
#include <Timers.d/GetTime.h>
#include <Corotational.d/GeomState.h>

#include <Control.d/ControlInterface.h>
#include <Solvers.d/Rbm.h>
#include <Utils.d/BinFileHandler.h>

#include <Utils.d/linkfc.h>
#include <Driver.d/GeoSource.h>

typedef FSFullMatrix FullMatrix;

template <class Scalar>
SingleDomainDynamic<Scalar>::SingleDomainDynamic(Domain *d)
{ 
  domain = d; 
  kelArray = 0; 
  allCorot = 0;
  geomState = 0; 
  userDefineDisp = 0;
  dprev = 0;

  flExchanger = domain->getFileExchanger();
}

//#define DEBUG_RBM_FILTER

template <class Scalar>
void
SingleDomainDynamic<Scalar>::projector_prep(Rbm *rbms, SparseMatrix *M)
{
 numR = rbms->numRBM();

 if (!numR) return;

 fprintf(stderr," ... Building the RBM Projector     ...\n");

 fprintf(stderr," ... Number of RBM(s)     =   %d     ...\n",numR);

 // KHP: store this pointer to the RBMs to use in the actual
 //      projection step within the time loop.

 int ndof = M->dim();
 Rmem = new double[numR*ndof];
 rbms->getRBMs(Rmem);
 StackFSFullMatrix Rt(numR, ndof, Rmem); 

 //DEBUG
//#ifdef DEBUG_RBM_FILTER
 FullMatrix R = Rt.transpose();
 //R.print("","R");
//#endif

 double *MRmem = new double[numR*ndof];
 StackFSFullMatrix MR(numR, ndof, MRmem);

 int n;
 for(n=0; n<numR; ++n)
   M->mult(Rmem+n*ndof, MRmem+n*ndof);

 FullMatrix MRt = MR.transpose();

 FullMatrix RtMR(numR,numR);
 Rt.mult(MRt,RtMR); 

 FullMatrix RtMRinverse = RtMR.invert();

 X = new FullMatrix(ndof,numR); 
 MRt.mult(RtMRinverse,(*X));

}

template <class Scalar>
void
SingleDomainDynamic<Scalar>::projector_prep(Rbm *rbms, GenSparseMatrix<DComplex> *M)
{
 numR = rbms->numRBM();
 if (!numR) return;
 fprintf(stderr," DynamDescr.C:projector_prep not implemented for complexM\n");
}

template <class Scalar>
void
SingleDomainDynamic<Scalar>::trProject(Vector &f)
{
 if (!numR) return;

 int ndof = f.size();

 double *yMem = (double *) dbg_alloca(numR*sizeof(double));
 double *zMem = (double *) dbg_alloca(ndof*sizeof(double));
 StackVector y(numR,yMem);
 StackVector z(ndof,zMem);

 StackFSFullMatrix Rt(numR, ndof, Rmem);

 // y = Rt*f
 Rt.mult(f,y);

 // z = X*y
 (*X).mult(y,z);

 // f = f - z;
 f.linC(1.0, f, -1.0, z);
}

template <class Scalar>
void
SingleDomainDynamic<Scalar>::project(Vector &v)
{
 if (!numR) return;

 int ndof = v.size();

 double *yMem = (double *) dbg_alloca(numR*sizeof(double));
 double *zMem = (double *) dbg_alloca(ndof*sizeof(double));
 StackVector y(numR,yMem);
 StackVector z(ndof,zMem);

 StackFSFullMatrix Rt(numR, ndof, Rmem);

 // y = Xt*v
 (*X).trMult(v,y);
 

 // z = R*y
 Rt.trMult(y,z);

 // v = v - z;
 v.linAdd(-1.0, z);
 (*X).trMult(v,y);
}

template <class Scalar>
int
SingleDomainDynamic<Scalar>::getFilterFlag()
{
 return domain->solInfo().filterFlags;
}

template <class Scalar>
int
SingleDomainDynamic<Scalar>::solVecInfo()
{
 return domain->numUncon();
}

template <class Scalar>
int
SingleDomainDynamic<Scalar>::dbcVecInfo()
{
 return domain->numdof();
}

template <class Scalar>
int
SingleDomainDynamic<Scalar>::bcInfo()
{
 return domain->nDirichlet();
}

template <class Scalar>
int
SingleDomainDynamic<Scalar>::getTimeIntegration()
{
 return domain->solInfo().timeIntegration;
}

template <class Scalar>
void
SingleDomainDynamic<Scalar>::getTimes(double &dt, double &tmax)
{
 dt   = domain->solInfo().dt;	// time step size
 tmax = domain->solInfo().tmax;	// total time
 if(userSupFunc)
   userSupFunc->setDt(dt);
}

template <class Scalar>
void
SingleDomainDynamic<Scalar>::getInitialTime(int &initTimeIndex, double &initTime)
{
 initTimeIndex = domain->solInfo().initialTimeIndex;
 initTime      = domain->solInfo().initialTime;
 t0 = initTime;
}

template <class Scalar>
double
SingleDomainDynamic<Scalar>::getInitialForceNorm()
{
  return domain->solInfo().initExtForceNorm;
}

template <class Scalar>
void
SingleDomainDynamic<Scalar>::getNewMarkParameters(double &beta, double &gamma,
                                          double &alphaf, double &alpham)
{
 beta  = domain->solInfo().newmarkBeta;
 gamma = domain->solInfo().newmarkGamma;
 alphaf = domain->solInfo().newmarkAlphaF;
 alpham = domain->solInfo().newmarkAlphaM;
}

template <class Scalar>
void
SingleDomainDynamic<Scalar>::getQuasiStaticParameters(double &maxVel, double &delta) 
{
 maxVel  = domain->solInfo().qsMaxvel;
 delta  = domain->solInfo().delta;
}

template <class Scalar>
void
SingleDomainDynamic<Scalar>::getRayleighCoef(double &alpha)
{
 alpha = domain->solInfo().alphaDamp;
}

template <class Scalar>
void
SingleDomainDynamic<Scalar>::getSteadyStateParam(int &steadyFlag, int &steadyMin, 
                                         int &steadyMax, double &steadyTol)
{
 steadyFlag  = domain->solInfo().steadyFlag;
 steadyMin   = domain->solInfo().steadyMin;
 steadyMax   = domain->solInfo().steadyMax;
 steadyTol   = domain->solInfo().steadyTol;
}

template <class Scalar>
void
SingleDomainDynamic<Scalar>::makeInvDiagMass(SparseMatrix *M, Vector& mdiag)
{
 int i;
 for(i=0; i<M->dim(); ++i)
   mdiag[i] = 1.0 / M->diag(i);
}

template <class Scalar>
void
SingleDomainDynamic<Scalar>::computeStabilityTimeStep(double& dt, DynamMat& dMat)
{
 // ... Compute Stability Time Step
 double sts = domain->computeStabilityTimeStep(dMat);

 fprintf(stderr," **************************************\n");
 if (domain->solInfo().modifiedWaveEquation) {
   sts = 1.73205*sts;
   fprintf(stderr," CONDITIONALLY STABLE MODIFIED WAVE EQUATION\n");
 }
 else
   fprintf(stderr," CONDITIONALLY STABLE NEWMARK ALGORITHM\n");
 fprintf(stderr," --------------------------------------\n");
 fprintf(stderr," Specified time step      = %10.4e\n",dt);
 fprintf(stderr," Stability max. time step = %10.4e\n",sts);
 fprintf(stderr," **************************************\n");
 if( sts < dt ) {
   dt = sts;
   fprintf(stderr," Stability max. time step is selected\n");
 } else
   fprintf(stderr," Specified time step is selected\n");
 fprintf(stderr," **************************************\n");
 fflush(stderr);

 domain->solInfo().setTimes(domain->solInfo().tmax, dt, domain->solInfo().dtemp);
}

template <class Scalar>
void
SingleDomainDynamic<Scalar>::getInitState(SysState<Vector> &inState)
{
  // initialize state with IDISP/IDISP6/IVEL/IACC or RESTART (XXXX initial accelerations are currently not supported)
  if(domain->solInfo().order == 1)
    domain->initTempVector(inState.getDisp(), inState.getVeloc(), inState.getPrevVeloc());
  else
    domain->initDispVeloc(inState.getDisp(),  inState.getVeloc(),
                          inState.getAccel(), inState.getPrevVeloc()); // IVEL, IDISP, IDISP6, restart

  // if we have a user supplied function, give it the initial state at the sensors
  // .. first update bcx, vcx in case any of the sensors have prescribed displacements
  if(claw && userSupFunc) {
    if(claw->numUserDisp) {
      double *userDefineDisp = new double[claw->numUserDisp];
      double *userDefineVel = new double[claw->numUserDisp];
      userSupFunc->usd_disp(domain->solInfo().initialTime, userDefineDisp, userDefineVel);
      setBC(userDefineDisp, userDefineVel);
      delete [] userDefineDisp; delete [] userDefineVel;
    }
    if(claw->numSensor) {
      double *ctrdisp = new double[claw->numSensor];
      double *ctrvel = new double[claw->numSensor];
      double *ctracc = new double[claw->numSensor];
      extractControlData(inState, ctrdisp, ctrvel, ctracc);
      userSupFunc->init(ctrdisp, ctrvel, ctracc, this);
      delete [] ctrdisp; delete [] ctrvel; delete [] ctracc;
    }
  }
}

template <class Scalar>
void
SingleDomainDynamic<Scalar>::getPitaState(SysState<Vector> &inState, int sliceRank)
{
  // PITA: Use the user-supplied initial seed values
  if (sliceRank == 0)
  {
    getInitState(inState);
  }
  else
  {
    Vector  &d_n = inState.getDisp();
    Vector  &v_n = inState.getVeloc();
    Vector  &v_p = inState.getPrevVeloc();
    Vector  &a_n = inState.getAccel();
    d_n.zero(); v_n.zero(); a_n.zero(); v_p.zero();
    domain->initDispVelocOnTimeSlice(d_n, v_n, sliceRank);

    // If we have a user supplied function, extract user supplied
    // initial displacement and velocity
    if (userSupFunc)
    {
      double sliceTime = domain->solInfo().dt * domain->solInfo().Jratio * sliceRank;
      double *ctrdisp = (double *) dbg_alloca(sizeof(double) * claw->numSensor);
      double *ctrvel  = (double *) dbg_alloca(sizeof(double) * claw->numSensor);
      double *ctracc  = (double *) dbg_alloca(sizeof(double) * claw->numSensor);
      extractControlData(inState, ctrdisp, ctrvel, ctracc);
      userSupFunc->usd_disp(sliceTime, ctrdisp, ctrvel);
    }
  }
}

template <class Scalar>
void
SingleDomainDynamic<Scalar>::extractControlData(SysState<Vector> &state, double *ctrdsp, 
                                                double *ctrvel, double *ctracc)
{
  // get SENSOR state
  Vector &dsp = state.getDisp();
  Vector &vlc = state.getVeloc();
  Vector &acc = state.getAccel();

  DofSetArray *cdsa = domain->getCDSA();
  DofSetArray *dsa = domain->getDSA();

  for(int i = 0; i < claw->numSensor; ++i) {
    int dof = cdsa->locate(claw->sensor[i].nnum, 1 << claw->sensor[i].dofnum);
    if(dof >= 0) { // free
      ctrdsp[i] = dsp[dof];
      ctrvel[i] = vlc[dof];
      ctracc[i] = acc[dof];
    }
    else { // either constrained or non-existant
      int dof2 = dsa->locate(claw->sensor[i].nnum, 1 << claw->sensor[i].dofnum);
      if(dof2 >= 0) { // constrained
        ctrdsp[i] = bcx[dof2];
        ctrvel[i] = vcx[dof2];
        ctracc[i] = 0.0; // XXXX prescribed acceleration not supported
      }
    }
  }
}

template <class Scalar>
void
SingleDomainDynamic<Scalar>::addCtrl(Vector &f, double *ctrfrc)
{
  // add ACTUATOR forces
  DofSetArray *cdsa = domain->getCDSA();
  for(int i = 0; i < claw->numActuator; ++i) {
    int dof = cdsa->locate(claw->actuator[i].nnum, 1 << claw->actuator[i].dofnum);
    if(dof >= 0) f[dof] += ctrfrc[i];
  }
}

template <class Scalar>
void
SingleDomainDynamic<Scalar>::addUserForce(Vector&f, double *userDefineForce)
{
  // add USDF forces
  DofSetArray *cdsa = domain->getCDSA();
  for(int i = 0; i < claw->numUserForce; ++i) {
    int dof = cdsa->locate(claw->userForce[i].nnum, 1<<claw->userForce[i].dofnum);
    if(dof >= 0) f[dof] += userDefineForce[i];
  }
}

template <class Scalar>
void
SingleDomainDynamic<Scalar>::setBC(double *userDefineDisp, double* userDefineVel)
{
  // update time-dependent prescribed displacements and velocities
  DofSetArray *dsa = domain->getDSA();
  for(int i = 0; i < claw->numUserDisp; ++i) {
    int dof = dsa->locate(claw->userDisp[i].nnum,1 << claw->userDisp[i].dofnum);
    if(dof >= 0) {
      bcx[dof] = userDefineDisp[i];
      vcx[dof] = userDefineVel[i];
    }
  }
}

template <class Scalar>
void
SingleDomainDynamic<Scalar>::getConstForce(Vector &cnst_f)
{
  cnst_f.zero(); // constant force i.e. gravity + non-follower pressure

  // gravity force computation
  if(domain->gravityFlag()) domain->template buildGravityForce<double>(cnst_f);

  // pressure force is also constant with respect to time for linear analysis (by convention)
  if(domain->pressureFlag() && !domain->solInfo().isNonLin()) domain->template buildPressureForce<double>(cnst_f);
}

template <class Scalar>
void
SingleDomainDynamic<Scalar>::getContactForce(Vector &d, Vector &ctc_f)
{
  times->tdenforceTime -= getTime();
  ctc_f.zero();
  if(domain->tdenforceFlag()) {
    times->updateSurfsTime -= getTime();
    domain->UpdateSurfaces(geomState, 1); // update to current configuration
    times->updateSurfsTime += getTime();

    Vector dinc(domain->numUncon());
    dinc.linC(1.0, d, -1.0, *dprev);
    geomState->update(dinc);
    times->updateSurfsTime -= getTime();
    domain->UpdateSurfaces(geomState, 2); // update to predicted configuration
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
void
SingleDomainDynamic<Scalar>::computeExtForce2(SysState<Vector> &state, Vector &ext_f, 
                                              Vector &cnst_f, int tIndex, double t,
                                               Vector *aero_f, double gamma, double alphaf, double *pt_dt)
{
  times->formRhs -= getTime();

  ext_f.zero();

  // compute USDD prescribed displacements
  double *userDefineDisp = 0;
  if(claw && userSupFunc) {
    if(claw->numUserDisp) {
      userDefineDisp = new double[claw->numUserDisp];
      double *userDefineVel = new double[claw->numUserDisp];
      userSupFunc->usd_disp(t, userDefineDisp, userDefineVel);
      setBC(userDefineDisp, userDefineVel); // update bcx, vcx
      delete [] userDefineVel;
    }
  }

  // update geomState for nonlinear problems. note computeExtForce2 must be called before getInternalForce
  if(domain->solInfo().isNonLin() || domain->tdenforceFlag()) {
    if(!dprev) { dprev = new Vector(domain->numUncon(), 0.0); }
    Vector dinc(domain->numUncon());
    dinc.linC(1.0, state.getDisp(), -1.0, *dprev); // incremental displacement: dinc = d - dprev
    geomState->update(dinc);
    if(userDefineDisp) geomState->updatePrescribedDisplacement(userDefineDisp, claw, domain->getNodes());
    (*dprev) = state.getDisp(); // keep a copy so the incremental displacement can be computed
    geomState->setVelocity(state.getVeloc()); 
  }

  // add FORCE (including MFTT), HDNB, ROBIN, GRAVITY and AEROELASTIC forces
  // for linear problems also add PRESSURE and non-homogeneous DISP/TEMP and USDD forces
  //int *nomap = 0;
  domain->computeExtForce4(*prevFrc, ext_f, cnst_f, tIndex, t, 
                           kuc, userDefineDisp, (int *) 0, aero_f, gamma, alphaf);
  if(userDefineDisp) delete [] userDefineDisp;

  // for nonlinear problems add PRESSURE forces
  if(domain->solInfo().isNonLin() && domain->pressureFlag()) { 
    domain->buildPressureForce(ext_f, geomState);
  }

  // add USDF forces
  if(claw && userSupFunc) {
    if(claw->numUserForce) {
      double *userDefineForce = new double[claw->numUserForce];
      userSupFunc->usd_forc(t, userDefineForce);
      addUserForce(ext_f, userDefineForce);
      delete [] userDefineForce;
    }
  }

  // add ACTUATOR forces XXXX but is state at time t? need to check disp/vel/acc for explicit/implicit/quasistatics
  if(claw && userSupFunc) {
    if(claw->numActuator > 0) {
      double *ctrdisp = new double[claw->numSensor];
      double *ctrvel = new double[claw->numSensor];
      double *ctracc = new double[claw->numSensor];
      double *ctrfrc = new double[claw->numActuator];

      extractControlData(state, ctrdisp, ctrvel, ctracc);
      userSupFunc->ctrl(ctrdisp, ctrvel, ctracc, ctrfrc, t, &state, &ext_f);
      addCtrl(ext_f, ctrfrc);
      delete [] ctrdisp; delete [] ctrvel; delete [] ctracc; delete [] ctrfrc;
    }
  }
 
  // KHP: apply projector here
  int useProjector = domain->solInfo().filterFlags;
  if(useProjector) trProject(ext_f); 

  if(tIndex == 1)
    domain->solInfo().initExtForceNorm = ext_f.norm();

  times->formRhs += getTime();
}

template <class Scalar>
void
SingleDomainDynamic<Scalar>::processLastOutput()  {

  OutputInfo *oinfo = geoSource->getOutputInfo();
  domain->solInfo().lastIt = true;
  for (int iOut = 0; iOut < geoSource->getNumOutInfo(); iOut++)
    oinfo[iOut].interval = 1;
}

template <class Scalar>
void
SingleDomainDynamic<Scalar>::preProcess()
{
  // Allocate space for Timers
  times = new StaticTimers;
  startTimerMemory(times->preProcess, times->memoryPreProcess);

  // Makes renumbering, connectivities and dofsets
  domain->preProcessing();

  // initialize bcx and vcx with time-independent prescribed displacements and velocities
  // note: bcx and vcx are updated with time-dependent prescribed displacements and velocities (USDD) later
  times->makeBCs -= getTime();
  int numdof = domain->numdof();
  int *bc = new int[numdof];
  bcx = new double[numdof];
  domain->make_bc(bc, bcx);
  delete [] bc;
  vcx = new double[numdof]; 
  for(int i=0; i<numdof; ++i) vcx[i] = 0.0;
  times->makeBCs += getTime();

  times->makeDOFs -= getTime();
  domain->make_constrainedDSA();
  domain->makeAllDOFs();
  times->makeDOFs += getTime();

  // Check for user supplied routines (control, force or displacement)
  claw = geoSource->getControlLaw();
  userSupFunc = domain->getUserSuppliedFunction();

  if((domain->numInitDisp6() > 0 && domain->solInfo().gepsFlg == 1) || domain->solInfo().isNonLin()) { // PJSA 3-31-08
    FullSquareMatrix *geomKelArray=0;
    // this function builds corotators, geomstate and kelArray 
    // for linear+geps only it updates geomState with IDISP6 and computes the element stiffness matrices using this updated geomState
    domain->computeGeometricPreStress(allCorot, geomState, kelArray, times, geomKelArray);
  }
  else if(domain->tdenforceFlag()) geomState = new GeomState(*domain->getDSA(), *domain->getCDSA(), domain->getNodes());

  if(domain->solInfo().isNonLin() || domain->tdenforceFlag()) {
    // for nonlinear explicit we only need to initialize geomState with the constant constrained displacements (DISP).
    // the geomState is always updated before use with the current unconstrained displacments plus any time-dependent constrained displacements (USDD)
    if(domain->nDirichlet() > 0) { 
      geomState->updatePrescribedDisplacement(domain->getDBC(), domain->nDirichlet()); 
    }
  }

  if(domain->tdenforceFlag()) {
    domain->InitializeDynamicContactSearch();
  }

  prevFrc = new PrevFrc(domain->numUncon());
  prevFrcBackup = new PrevFrc(domain->numUncon());

  stopTimerMemory(times->preProcess, times->memoryPreProcess);
}

// In general a symmetric matrix can be partitioned as:
//
// A = | A11   A12 |
//     | A12^t A22 |
//

template <class Scalar>
GenDynamMat<Scalar> *
SingleDomainDynamic<Scalar>::buildOps(double coeM, double coeC, double coeK)
{

 // KHP: 7-30-98 added sparse matrices for Muc and Cuc 
 // to allow for prescribed displacements
 // that are changing with respect to time.

 AllOps<Scalar> allOps;
 GenDynamMat<Scalar> *dMat = new GenDynamMat<Scalar>;

 allOps.K   = domain->template constructDBSparseMatrix<Scalar>();
 allOps.M   = domain->template constructDBSparseMatrix<Scalar>();
 allOps.Muc = domain->template constructCuCSparse<Scalar>();
 allOps.Kuc = domain->template constructCuCSparse<Scalar>();
 allOps.Mcc = domain->template constructCCSparse<Scalar>();

 // Rayleigh Damping coefficients
 double alpha = domain->solInfo().alphaDamp;
 double beta  = domain->solInfo().betaDamp;

 // Damping Matrix: C = alpha*M + beta*K + D
 // note #1: for explicit central difference time integration (newmarkBeta = 0.0) rayleigh mass damping is embedded
 //          and therefore it is not necessary to assemble the C matrix.
 if((alpha != 0.0 && domain->solInfo().newmarkBeta != 0.0) || beta != 0.0 || domain->solInfo().ATDARBFlag != -2.0) {
   allOps.C   = domain->template constructDBSparseMatrix<Scalar>();
   allOps.Cuc = domain->template constructCuCSparse<Scalar>();
 }

 // to compute a^0 = M^{-1}(f_ext^0-f_int^0-Cu^0)
 if(domain->solInfo().newmarkBeta != 0.0 && domain->solInfo().iacc_switch) { // not required for explicit
   GenBLKSparseMatrix<Scalar> *m = domain->template constructBLKSparseMatrix<Scalar>(domain->getCDSA());
   m->zeroAll();
   allOps.Msolver = m;
   dMat->Msolver = m;
 }

 Rbm *rigidBodyModes = 0;

 int useProjector = domain->solInfo().filterFlags;
 int useGrbm = domain->solInfo().rbmflg; 

 if (useGrbm || useProjector) 
   rigidBodyModes = domain->constructRbm();
   
 if (!useGrbm || (getTimeIntegration() != 1) ) 
   domain->template buildOps<Scalar>(allOps, coeK, coeM, coeC, 0, kelArray);
 else
   domain->template buildOps<Scalar>(allOps, coeK, coeM, coeC, rigidBodyModes, kelArray); 

 // Filter RBM
 if(useProjector == 1)
    fprintf(stderr," ... RBMfilter Level 1 Requested    ...\n");
 if(useProjector == 2)
    fprintf(stderr," ... RBMfilter Level 2 Requested    ...\n");

 if(useProjector)
   projector_prep(rigidBodyModes, allOps.M);
 
 // Modal decomposition preprocessing
 int decompFlag = domain->solInfo().modeDecompFlag;
 if(decompFlag) {
   fprintf(stderr," ... Modal decomposition requested ...\n");
   modeDecompPreProcess(allOps.M);
 }

 dMat->K      = allOps.K;
 dMat->M      = allOps.M;
 dMat->C      = allOps.C;
 dMat->Cuc    = allOps.Cuc;
 dMat->Muc    = allOps.Muc;
 dMat->Mcc    = allOps.Mcc;
 kuc          = allOps.Kuc;
 dMat->kuc    = allOps.Kuc;
 dMat->dynMat = allOps.sysSolver;
 if(dMat->Msolver) dMat->Msolver->factor();

 if(domain->tdenforceFlag()) domain->MakeNodalMass(allOps.M); 

 return dMat;
}

template <class Scalar>
void
SingleDomainDynamic<Scalar>::getInternalForce(Vector& d, Vector& f, double t)
{
  if(domain->solInfo().isNonLin()) { // PJSA 3-31-08
    Vector residual(domain->numUncon(),0.0);
    Vector fele(domain->maxNumDOF());
    domain->getStiffAndForce(*geomState, fele, allCorot, kelArray, residual); // kelArray = tangent stiffness, residual -= internal force (nonlinear elements)
    f.linC(-1.0,residual); // f = -residual
  }
  else {
    f.zero();
    domain->getKtimesU(d, bcx, f, 1.0, kelArray);  // note: although passed as an argument, the bcx contribution is not computed in this function
  }
  int useProjector = domain->solInfo().filterFlags;
  if(useProjector) trProject(f);
}

template <class Scalar>
void
SingleDomainDynamic<Scalar>::computeTimeInfo()
{
  // Time integration information

  // Get total time and time step size and store them
  double totalTime = domain->solInfo().tmax;
  double dt        = domain->solInfo().dt;

  // Compute maximum number of steps
  int maxStep = (int) ( (totalTime+0.49*dt)/dt );

  // Compute time remainder
  double remainder = totalTime - maxStep*dt;
  if(std::abs(remainder)>0.01*dt){
    domain->solInfo().tmax = maxStep*dt;
    fprintf(stderr, " Warning: Total time is being changed to : %e\n", domain->solInfo().tmax);
  }
}

template <class Scalar>
int 
SingleDomainDynamic<Scalar>::aeroPreProcess(Vector& d_n, Vector& v_n, 
                                            Vector& a_n, Vector& v_p)
{
  domain->aeroPreProcess(d_n, v_n, a_n, v_p, bcx, vcx);
  return domain->solInfo().aeroFlag;
}

template <class Scalar>
void
SingleDomainDynamic<Scalar>::a5TimeLoopCheck(int& parity, double& t, double dt)
{
  if(domain->solInfo().aeroFlag == 5) {
     if(!parity) t -= dt;
     parity = ( parity ? 0 : 1 );
  }
}

template <class Scalar>
void
SingleDomainDynamic<Scalar>::a5StatusRevise(int parity, SysState<Vector>& curState, 
                                            SysState<Vector>& bkState)
{
  if(domain->solInfo().aeroFlag == 5) {
     if(parity) { // restore

       *prevFrc = *prevFrcBackup;
       curState.getDisp()      = bkState.getDisp();
       curState.getVeloc()     = bkState.getVeloc();
       curState.getAccel()     = bkState.getAccel();
       curState.getPrevVeloc() = bkState.getPrevVeloc();

     } else {      // backup

       *prevFrcBackup = *prevFrc;
       bkState.getDisp()      = curState.getDisp();
       bkState.getVeloc()     = curState.getVeloc();
       bkState.getAccel()     = curState.getAccel();
       bkState.getPrevVeloc() = curState.getPrevVeloc();

     }
  }
}

template <class Scalar>
void
SingleDomainDynamic<Scalar>::thermoePreProcess(Vector&, Vector&, Vector&)
{
  domain->thermoePreProcess();
}


// Single Domain Post Processor Functions

template <class Scalar>
SDDynamPostProcessor *
SingleDomainDynamic<Scalar>::getPostProcessor()
{
  return new SDDynamPostProcessor(domain,bcx,vcx,times);
}

template <class Scalar>
void
SingleDomainDynamic<Scalar>::printTimers(DynamMat *dynamMat)
{

 long memoryUsed = (dynamMat->dynMat)->size();//solver->size();
 double solveTime  = dynamMat->dynMat->getSolutionTime();//solver->getSolutionTime();

 times->printStaticTimers(solveTime, memoryUsed, domain);

  if(domain->solInfo().massFlag)  {
   double mass = domain->computeStructureMass();
   fprintf(stderr," --------------------------------------\n");
   fprintf(stderr," ... Structure mass = %e  ...\n",mass);
   fprintf(stderr," --------------------------------------\n");
 }

}

template <class Scalar>
void
SingleDomainDynamic<Scalar>::addPrescContrib(SparseMatrix *Muc, SparseMatrix *Cuc,
                                     Vector& dnc, Vector& vnc, Vector& anc,
                                     Vector& result, double t, double *pt_dt)
{

 double dt    = domain->solInfo().dt;
 double beta  = domain->solInfo().newmarkBeta;
 double gamma = domain->solInfo().newmarkGamma;

 // Compute acceleration at constrained degrees of freedom
 anc *= (1.0 - 1.0/(2.0*beta));
 anc.linAdd( -(1.0/(dt*dt*beta)), dnc, -(1.0/(dt*beta)), vnc );

 // prescribed displacement boundary conditions at half time step
 Vector dn_h(dnc);

 double *userDefineDisp = 0;
 double *userDefineVel  = 0;

 int i;
 if( claw && userSupFunc ) {

   userDefineDisp = (double *) dbg_alloca( sizeof(double)*claw->numUserDisp );

   userDefineVel  = (double *) dbg_alloca( sizeof(double)*claw->numUserDisp );

   userSupFunc->usd_disp( t+dt/2.0, userDefineDisp, userDefineVel );

   for(i=0; i<claw->numUserDisp; ++i) {
     int dof = domain->getDSA()->locate( claw->userDisp[i].nnum,
                                         1 << claw->userDisp[i].dofnum );
     if(dof < 0) continue;
     int dof1 = domain->getCDSA()->invRCN( dof );
     if(dof1 >= 0)
       dn_h[dof1] = userDefineDisp[i];
   }

   userSupFunc->usd_disp( t, userDefineDisp, userDefineVel );

   for(i=0; i<claw->numUserDisp; ++i) {
     int dof = domain->getDSA()->locate( claw->userDisp[i].nnum,
                                         1 << claw->userDisp[i].dofnum );
     if(dof < 0) continue;
     int dof1 = domain->getCDSA()->invRCN( dof );
     if(dof1 >= 0) {
       dnc[dof1] = userDefineDisp[i];
       vnc[dof1] = userDefineVel[i];
     }
   }
  
 }

 // update acceleration at constrained points
 anc.linAdd( 1.0/(dt*dt*beta), dnc );

 Vector dis( dnc );

 // dis = dnc + 0.5*dt*vnc + dt*dt*(0.25 - beta)*anc
 dis.linAdd( 0.5*dt, vnc, dt*dt*(0.25 - beta), anc);
 
 Muc->mult( dis, result );

 // Viscous Damping Matrix contributions

 if( Cuc ) {

   // dis = (dt*gamma) * dnc - (dt*dt*(beta - 0.5*gamma)) * vnc 
   //         - (dt*dt*dt*(0.5*beta - 0.25*gamma)) * anc
   //         - (dt*gamma) * dn_h

   dis.linC( dt*gamma, dnc, -(dt*dt*(beta - 0.5*gamma)), vnc );

   dis.linAdd( -(dt*dt*dt*(0.5*beta - 0.25*gamma)), anc );

   dis.linAdd( -dt*gamma, dn_h );

   dis *= -1.0;

   Cuc->mult( dis, result );

 }

}

template <class Scalar>
double
SingleDomainDynamic<Scalar>::betaDamp() const {
  return domain->solInfo().betaDamp;
}

template <class Scalar>
double
SingleDomainDynamic<Scalar>::alphaDamp() const {
  return domain->solInfo().alphaDamp;
}

template <class Scalar>
void 
SingleDomainDynamic<Scalar>::setDamping( double betaDamp, double alphaDamp )
{ 
  domain->solInfo().setDamping( betaDamp, alphaDamp ); 
}

template <class Scalar>
int
SingleDomainDynamic<Scalar>::getModeDecompFlag()
{
 return domain->solInfo().modeDecompFlag;
}   

template <class Scalar>
void
SingleDomainDynamic<Scalar>::modeDecompPreProcess(SparseMatrix *M)
{

// Read Eigenvectors from file EIGENMODES
// ======================================
    int eigsize;
 
    BinFileHandler modefile("EIGENMODES" ,"r");
    modefile.read(&maxmode, 1);
    fprintf(stderr,"Number of Modes = %d\n", maxmode);
 
    modefile.read(&eigsize, 1);
    fprintf(stderr,"Size of EigenVector = %d\n", eigsize);
 
    eigmodes = new double*[maxmode];
    int i;
    for (i=0; i<maxmode; ++i)
     eigmodes[i] = new double[eigsize];

// Check if the problem sizes are identical

    int numdof = solVecInfo();

    if (eigsize != numdof)
      fprintf(stderr, "... ERROR !! EigenVector and Problem sizes differ... ");
      fflush(stderr);

    for (i=0; i<maxmode; ++i)
      modefile.read(eigmodes[i], eigsize);  

// Compute Transpose(Phi_i)*M once and for all
// ===========================================
 
   fprintf(stderr, "Preparing Transpose(Phi_i)*M\n");
 
   tPhiM =  new double*[maxmode];
     for (i=0; i<maxmode; ++i)
       tPhiM[i] = new double[eigsize];
 
   for (i = 0 ; i<maxmode; ++i)
     M->mult(eigmodes[i], tPhiM[i]);  // taking advantage of symmetry of M and computing
                                      //   M*Phi_i instead of transpose(Phi_i)*M
/*  DEBUG: UNCOMMENT FOR DEBUGGING

// Verify that Phi_i*M*Phi_i = 1 and Phi_i*M*Phi_j = 0
// ===================================================

   double PhiMPhi = 0;

   for (i =0; i<maxmode; ++i) {
     for (int j =0; j<maxmode; ++j) {
      for (int k = 0; k < eigsize; ++k) {
          PhiMPhi = PhiMPhi + eigmodes[i][k]*M->diag(k)*eigmodes[j][k];
      }
     fprintf(stderr, "Phi_%d*M*Phi_%d = %19.11e\n",i,j, PhiMPhi) ;
     PhiMPhi = 0;
     }
     fprintf(stderr,"\n");
   }

*/

//  Need eigmodes for error computation
//   delete eigmodes;

}

template <class Scalar>
void
SingleDomainDynamic<Scalar>::modeDecompPreProcess(GenSparseMatrix<DComplex> *M)
{
  fprintf(stderr,"DynamDesc.C::modeDecompPreProcess not implemented for complexM\n");
}

template <class Scalar>
void
SingleDomainDynamic<Scalar>::modeDecomp(double t, int tIndex, Vector& d_n)
{
// Compute Alpha and error only if their output file is specified,
// otherwise, it wouldn't make sense

   // Compute alfa_i=PhiDiag_i*d_n
 
  int numOutInfo = geoSource->getNumOutInfo();
  OutputInfo *oinfo = geoSource->getOutputInfo();

  int i, j, k;
  alfa = 0;

  for (i=0; i < numOutInfo; i++) {
    if(oinfo[i].interval != 0 && tIndex % oinfo[i].interval == 0) {

      int w = oinfo[i].width;
      int p = oinfo[i].precision;

      switch(oinfo[i].type) {
        case OutputInfo::ModeAlpha: {
          //fprintf(stderr, "Computing alfa_i\n");

          if(!alfa) {
            alfa = new double[maxmode];

            for (k=0; k<maxmode; ++k)
              alfa[k] = 0;

            for (j = 0; j < maxmode; ++j)
              for (k = 0; k < d_n.size(); ++k)
                 alfa[j] += tPhiM[j][k]*d_n[k];
          }

          // Write alfa
          fprintf(oinfo[i].filptr, "%e  ", t);
          for(j=0; j<maxmode; ++j)
            fprintf(oinfo[i].filptr, "% *.*E ", w, p, alfa[j]);
            fprintf(oinfo[i].filptr, "\n");
     
          fflush(oinfo[i].filptr);
        } break;

        case OutputInfo::ModeError: {
          //fprintf(stderr, "Computing relative error\n");

          if(!alfa) {
            alfa = new double[maxmode];
  
            for (k=0; k<maxmode; ++k)
              alfa[k] = 0;
  
            for (j = 0; j < maxmode; ++j)
              for (k = 0; k < d_n.size(); ++k)
                alfa[j] += tPhiM[j][k]*d_n[k];
          }

          double sumerror = 0;
          double normerror = 0;
          double sumdisp = 0;
          double normdisp = 0;

          int ersize = d_n.size();

          double *sumalfa = new double[ersize];
          double *error   = new double[ersize];

          for (k=0; k < ersize; ++k) sumalfa[k] = 0;

          for (k=0; k < ersize; ++k)
            for (j=0; j < maxmode; ++j)
              sumalfa[k] += alfa[j]*eigmodes[j][k];


          for (j=0; j < ersize; ++j) {
            error[j] = d_n[j]-sumalfa[j];
            sumerror += error[j]*error[j];
            sumdisp += d_n[j]*d_n[j];
          }

          normdisp = sqrt(sumdisp);
          if (normdisp == 0.0)  normerror = 0.0;
          else normerror = sqrt(sumerror)/normdisp;

          // Write error
          fprintf(oinfo[i].filptr, "%e % *.*E\n", t, w, p, normerror);
          fflush(oinfo[i].filptr);
        } break;
   
        default: 
          break;
      }

    }
  }
  if(alfa) delete [] alfa;

}

template <class Scalar>
void
SingleDomainDynamic<Scalar>::aeroSend(double time, Vector& d, Vector& v, Vector& a, Vector& v_p)
{
  // XXXX need to compute bcx properly for USDD so appropriate prescribed displacements are sent to fluid
  domain->aeroSend(d, v, a, v_p, bcx, vcx);
}

