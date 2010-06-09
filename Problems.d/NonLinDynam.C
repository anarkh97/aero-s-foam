#include <Utils.d/dbg_alloca.h>
#include <stdio.h>
#include <stdlib.h>

#include <Driver.d/Domain.h>
#include <Problems.d/NonLinDynam.h>
#include <Problems.d/DynamDescr.h>
#include <Solvers.d/Solver.h>
#include <Timers.d/StaticTimers.h>
#include <Math.d/SparseMatrix.h>
#include <Math.d/CuCSparse.h>
#include <Math.d/DBSparseMatrix.h>
#include <Math.d/FullSquareMatrix.h>
#include <Control.d/ControlInterface.h>
#include <Corotational.d/TemperatureState.h>
#include <Corotational.d/GeomState.h>
#include <Corotational.d/Corotator.h>
#include <Timers.d/GetTime.h>
#include <Math.d/FullMatrix.h>
#include <Solvers.d/Rbm.h>
#include <Driver.d/GeoSource.h>
#include <Element.d/State.h>
#include <Driver.d/SysState.h>

typedef FSFullMatrix FullMatrix;

extern int verboseFlag;

NonLinDynamic::NonLinDynamic(Domain *d)
: domain(d),
  res((FILE*) 0),
  clawDofs(0),
  secondRes(0.0),
  numSystems(0)
{
  claw = 0; userSupFunc = 0;
}

NonLinDynamic::~NonLinDynamic()
{
  if (res != (FILE*) 0)
    fclose(res);
}

void
NonLinDynamic::projector_prep(Rbm *rbms, SparseMatrix *M)
{
 
 numR = rbms->numRBM();

 if (!numR) return;

 fprintf(stderr," ... Building the RBM Projector     ...\n");
 //fprintf(stderr," ... Number of RBM(s)     =   %d     ...\n",numR);

 int ndof = M->dim();

 // KHP: store this pointer to the RBMs to use in the actual
 //      projection step within the time loop.
 Rmem = new double[numR*ndof];
 rbms->getRBMs(Rmem);
 StackFSFullMatrix Rt(numR, ndof, Rmem);

 //DEBUG
#ifdef DEBUG_RBM_FILTER
 FullMatrix R = Rt.transpose();
 R.print("","R");
#endif

 double *MRmem = new double[numR*ndof];
 StackFSFullMatrix MR(numR, ndof, MRmem);

 int n;
 for(n=0; n<numR; ++n)
   M->mult(Rmem+n*ndof, MRmem+n*ndof);

 FullMatrix MRt = MR.transpose();
 //MRt.print("MR","MR");

 FullMatrix RtMR(numR,numR);
 Rt.mult(MRt,RtMR);
 //DEBUG
 //RtMR.print("RtMR","RtMR");

 FullMatrix RtMRinverse = RtMR.invert();
 //DEBUG
 //RtMRinverse.print("RtMRinverse");

 X = new FullMatrix(ndof,numR);
 MRt.mult(RtMRinverse,(*X));
 //DEBUG
 //X->print("MRRtMRinverse");

}

void
NonLinDynamic::trProject(Vector &f)
{
 if (!numR) return;

 int ndof = f.size();
 
 double *yMem = (double *) dbg_alloca(numR*sizeof(double));
 double *zMem = (double *) dbg_alloca(ndof*sizeof(double));
 StackVector y(ndof,yMem);
 StackVector z(ndof,zMem);

 StackFSFullMatrix Rt(numR, ndof, Rmem);

 // y = Rt*f
 Rt.mult(f,y);


 // z = X*y
 (*X).mult(y,z);


 // f = f - z;
 f.linC(1.0, f, -1.0, z);

}

void
NonLinDynamic::readRestartFile(Vector &d_n, Vector &v_n, Vector &a_n,
                               Vector &v_p, GeomState &geomState)
{
 domain->readRestartFile(d_n, v_n, a_n, v_p, bcx, vcx, geomState);
}

int
NonLinDynamic::getInitState(Vector& d_n, Vector& v_n, Vector &a_n, Vector &v_p)
{
  // initialize state with IDISP/IDISP6/IVEL/IACC or RESTART (XXXX initial accelerations are currently not supported)
  domain->initDispVeloc(d_n, v_n, a_n, v_p);

  updateUserSuppliedFunction(d_n, v_n, a_n, v_p, domain->solInfo().initialTime);

  int aeroAlg = domain->solInfo().aeroFlag; // by default, non-aeroelastic computation
  // call aeroPreProcess if a restart file does not exist
  if(aeroAlg >= 0 && geoSource->getCheckFileInfo()->lastRestartFile == 0) 
    domain->aeroPreProcess(d_n, v_n, a_n, v_p, bcx, vcx);

  if(domain->solInfo().aeroheatFlag >= 0)
    domain->aeroHeatPreProcess(d_n, v_n, v_p, bcx);

  if(domain->solInfo().thermoeFlag >= 0)
    domain->thermoePreProcess();

  if(domain->solInfo().thermohFlag >= 0)
    domain->thermohPreProcess(d_n, v_n, v_p, bcx);

  return aeroAlg;
}

void
NonLinDynamic::updateUserSuppliedFunction(Vector& d_n, Vector& v_n, Vector &a_n, Vector &v_p, double initialTime)
{
  // if we have a user supplied function, give it the initial state at the sensors
  // .. first update bcx, vcx in case any of the sensors have prescribed displacements
  if(claw && userSupFunc) {
    if(claw->numUserDisp) {
      double *userDefineDisp = new double[claw->numUserDisp];
      double *userDefineVel = new double[claw->numUserDisp];
      userSupFunc->usd_disp(initialTime, userDefineDisp, userDefineVel);
      setBC(userDefineDisp, userDefineVel);
      delete [] userDefineDisp; delete [] userDefineVel;
    }
    if(claw->numSensor) {
      double *ctrdisp = new double[claw->numSensor];
      double *ctrvel = new double[claw->numSensor];
      double *ctracc = new double[claw->numSensor];
      extractControlData(d_n, v_n, a_n, ctrdisp, ctrvel, ctracc); 
      userSupFunc->init(ctrdisp, ctrvel, ctracc);
      delete [] ctrdisp; delete [] ctrvel; delete [] ctracc;
    }
  }
}

void
NonLinDynamic::extractControlData(Vector& d_n, Vector& v_n, Vector& a_n,
                                  double *ctrdsp, double *ctrvel, double *ctracc)
{
  // get SENSOR states
  DofSetArray *cdsa = domain->getCDSA();
  DofSetArray *dsa = domain->getDSA();

  for(int i = 0; i < claw->numSensor; ++i) {
    int dof = cdsa->locate(claw->sensor[i].nnum, 1 << claw->sensor[i].dofnum);
    if(dof >= 0) { // free
      ctrdsp[i] = d_n[dof];
      ctrvel[i] = v_n[dof];
      ctracc[i] = a_n[dof];
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

void
NonLinDynamic::extractControlDisp(GeomState *geomState, double *ctrdsp)  
{
  CoordSet &nodes = domain->getNodes();
  NodeState *nodeState = geomState->getNodeState();
  for(int i = 0; i < claw->numSensor; ++i) {
    switch(claw->sensor[i].dofnum) {
      case 0:
        ctrdsp[i] = nodeState[claw->sensor[i].nnum].x - nodes[claw->sensor[i].nnum]->x;
        break;
      case 1:
        ctrdsp[i] = nodeState[claw->sensor[i].nnum].y - nodes[claw->sensor[i].nnum]->y;
        break;
      case 2:
        ctrdsp[i] = nodeState[claw->sensor[i].nnum].z - nodes[claw->sensor[i].nnum]->z;
        break;
      default:
        fprintf(stderr, "ERROR: Sensor dof %d not available in NonLinDynamic::extractControlDisp\n",claw->sensor[i].dofnum+1);
    }
  }
}

/*
void
NonLinDynamic::addCtrl(Vector& f, double *ctrfrc)
{
  DofSetArray *cdsa = domain->getCDSA();
  for(int i = 0; i < claw->numActuator; ++i) {
    int dof = cdsa->locate(claw->actuator[i].nnum,1 << claw->actuator[i].dofnum);
    if(dof >= 0) f[dof] += ctrfrc[i];
  }
}

void
NonLinDynamic::addUserForce(Vector& f, double *userDefinedForce)
{
  DofSetArray *cdsa = domain->getCDSA();
  for(int i = 0; i < claw->numUserForce; ++i) {
    int dof = cdsa->locate(claw->userForce[i].nnum,1<<claw->userForce[i].dofnum);
    if(dof >= 0) f[dof] += userDefinedForce[i];
  }
}
*/

void
NonLinDynamic::getConstForce(Vector& constantForce)
{
  domain->computeConstantForce(constantForce);
}

void
NonLinDynamic::computeTimeInfo()
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

  // set half time step size in user defined functions 
  if(userSupFunc)
    userSupFunc->setDt(delta);
}

double
NonLinDynamic::getStiffAndForce(GeomState& geomState, Vector& residual,
                                Vector& elementInternalForce, double t)
{
  times->buildStiffAndForce -= getTime();

  // update the geomState according to the USDD prescribed displacements
  if(claw && userSupFunc) {
    if(claw->numUserDisp > 0) {
      double *userDefineDisp = new double[claw->numUserDisp];
      double *userDefineVel  = new double[claw->numUserDisp];

      userSupFunc->usd_disp(t, userDefineDisp, userDefineVel); // XXXX should we do something with the userDefineVel?

      geomState.updatePrescribedDisplacement(userDefineDisp, claw, domain->getNodes());

      setBC(userDefineDisp, userDefineVel);
      delete [] userDefineDisp; delete [] userDefineVel;
    }

    if(claw->numActuator > 0) {
      double *ctrdisp = new double[claw->numSensor];
      double *ctrvel  = new double[claw->numSensor];
      double *ctracc  = new double[claw->numSensor];
      double *ctrfrc  = new double[claw->numActuator];

      for(int j = 0; j < claw->numSensor; j++) ctrvel[j] = ctracc[j] = 0.0; // XXXX f(v,a) currently not supported

      // KHP: we need the state of the control sensors to pass to
      //      the user supplied control function
      extractControlDisp(&geomState, ctrdisp);

      userSupFunc->ctrl(ctrdisp, ctrvel, ctracc, ctrfrc, t);
      domain->updateActuatorsInNbc(ctrfrc);

      delete [] ctrdisp; delete [] ctrvel; delete [] ctracc; delete [] ctrfrc;
    }
  }

  domain->getStiffAndForce(geomState, elementInternalForce, allCorot, kelArray, residual);

  times->buildStiffAndForce +=  getTime();
 
  // return residual force norm
  return residual.norm();
}

int
NonLinDynamic::getNumStages()
{
/*
 const double targetFirstResidual= 1.0E5;
 int numStages=2;
 if(firstRes > targetFirstResidual)
   numStages= int(firstRes/targetFirstResidual)+2;
 return numStages;
*/
 return int(0.2+ domain->solInfo().getNLInfo().maxLambda / domain->solInfo().getNLInfo().dlambda);
}

int
NonLinDynamic::checkConvergence(int iteration, double normRes, Vector &residual, Vector& dv, 
                                double time)
{
     /*if(dofTypeArray == 0)
       dofTypeArray = cdsa->makeDofTypeArray();*/
#ifdef PRINT_FORCENORMS
     ConstrainedDSA *cdsa = domain->getCDSA();
     double momenNorm = 0.0;
     double forceNorm = 0.0;
     int i;
     for(i=0; i<domain->numNodes(); ++i) {
       if(domain->solInfo().order == 1) {
         int dof = cdsa->locate(i, DofSet::Temp);
         if(dof >= 0)
           forceNorm += (residual[dof]*residual[dof]);
       }
       else {
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
     }
     if(iteration == 0) {
       firstForceNorm = (forceNorm == 0.0) ? 1.0 : sqrt(forceNorm);
       firstMomenNorm = (momenNorm == 0.0) ? 1.0 : sqrt(momenNorm);
     }
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
     if(iteration == 1) secondRes = normRes;

     double relRes = normRes/firstRes;
     double relDv  = normDv /firstDv;
     double relEng = normEnergy/firstEng;

     int converged = 0;

     if(normRes <= tolerance*firstRes) 
       converged = 1;

     // Check for divergence
     if(normRes >= 1.0e10 * firstRes && normRes > secondRes) {
       converged = -1;
     }

     if(verboseFlag) {
       fprintf(stderr," Iteration # %d\n",iteration);
       fprintf(stderr," r      = %e dv      = %e energy      = %e\n"
                      " rel. r = %e rel. dv = %e rel. energy = %e\n",
                        normRes,normDv,normEnergy,
                        relRes,relDv,relEng);
     }

     totIter++;
     fprintf(res,"%d %19.12e %e %e %e %e\n",totIter,time,normRes,relRes, normDv, relDv);
     fflush(res);

     // Store residual norm and dv norm for output
     times->norms[numSystems].normDv      = normDv;
     times->norms[numSystems].relativeDv  = relDv;
     times->norms[numSystems].normRes     = normRes;
     times->norms[numSystems].relativeRes = relRes;
     times->numSystems = numSystems;

     numSystems += 1;

     return converged;
}

GeomState*
NonLinDynamic::createGeomState()
{
  if(domain->solInfo().soltyp == 2)
    return new TemperatureState( *domain->getDSA(), *domain->getCDSA(), domain->getNodes() );
  else
    return new GeomState( *domain->getDSA(), *domain->getCDSA(), domain->getNodes() );
}

GeomState*
NonLinDynamic::copyGeomState(GeomState* geomState)
{
  if(domain->solInfo().soltyp == 2)
    return new TemperatureState(* (TemperatureState *) geomState);
  else
    return new GeomState(*geomState);
}

// Rebuild dynamic mass matrix
void
NonLinDynamic::reBuild(GeomState& geomState, int iteration, double localDelta)
{
 // note: localDelta = deltat/2
 times->rebuild -= getTime();

 // Rebuild every updateK iterations
 if(iteration % domain->solInfo().getNLInfo().updateK == 0)  {
   //fprintf(stderr, "Rebuilding Tangent Stiffness for Iteration %d\n", iteration);

   //PJSA 11/5/09: new way to rebuild solver (including preconditioner) and Kuc, now works for any solver
   spm->zeroAll();
   AllOps<double> ops;
   ops.Kuc = kuc;
   if(spp) {
     spp->zeroAll();
     ops.spp = spp;
   }
   double beta, gamma, alphaf, alpham, dt = 2*localDelta;
   getNewmarkParameters(beta, gamma, alphaf, alpham);
   double Kcoef = (domain->solInfo().order == 1) ? dt*gamma : dt*dt*beta;
   double Ccoef = (domain->solInfo().order == 1) ? 0 : dt*gamma;
   double Mcoef = (domain->solInfo().order == 1) ? 1 : (1-alpham)/(1-alphaf);
   domain->makeSparseOps<double>(ops, Kcoef, Mcoef, Ccoef, spm, kelArray, melArray);
   if(!verboseFlag) solver->setPrintNullity(false);
   solver->factor();
   if(prec) prec->factor();

 }
 times->rebuild += getTime();
}

int
NonLinDynamic::solVecInfo()
{
  return domain->numUncon();
}

int
NonLinDynamic::sysVecInfo()
{
  return domain->numdof();
}

int
NonLinDynamic::elemVecInfo()
{
 return domain->maxNumDOF();
}

int
NonLinDynamic::getMaxit()
{
 return domain->solInfo().getNLInfo().maxiter;
}

double
NonLinDynamic::getDeltaLambda()
{
 return domain->solInfo().getNLInfo().dlambda;
}

void
NonLinDynamic::getExternalForce(Vector& rhs, Vector& constantForce, int tIndex, double t, 
                                GeomState* geomState, Vector& elemNonConForce, 
                                Vector &aeroForce)
{
  // ... BUILD THE EXTERNAL FORCE at t_{n+1-alphaf}
  times->formRhs -= getTime();

  // update USDF
  if(claw && userSupFunc) {
    if(claw->numUserForce > 0) {
      double *userDefinedForce = new double[claw->numUserForce];
      userSupFunc->usd_forc(t, userDefinedForce);
      domain->updateUsdfInNbc(userDefinedForce);
      delete [] userDefinedForce;
    }
  }

  // update THERMOE 
  if(domain->solInfo().thermoeFlag >= 0 && tIndex >= 0)
    domain->thermoeComm();

  // add f(t) to constantForce (not including follower forces)
  domain->computeExtForce4(rhs, constantForce, t);

  // add aeroelastic forces from fluid dynamics code
  double beta, gamma, alphaf, alpham;
  getNewmarkParameters(beta, gamma, alphaf, alpham);
  if(domain->solInfo().aeroFlag >= 0 && tIndex >= 0) {
    domain->buildAeroelasticForce(rhs, *prevFrc, tIndex, t, gamma, alphaf);
  }

  // add aerothermal fluxes from fluid dynamics code
  if(domain->solInfo().aeroheatFlag >= 0 && tIndex >= 0)
    domain->buildAeroheatFlux(rhs, prevFrc->lastFluidLoad, tIndex, t);

 //  HAI: apply projector here
 if(domain->solInfo().filterFlags || domain->solInfo().hzemFilterFlag)
   trProject(rhs);

 times->formRhs += getTime();
}

void
NonLinDynamic::formRHSinitializer(Vector &fext, Vector &velocity, Vector &elementInternalForce, GeomState &geomState, Vector &rhs)
{
  // rhs = (fext - fint - Cv)
  rhs = fext;
  elementInternalForce.zero();
  domain->getStiffAndForce(geomState, elementInternalForce, allCorot, kelArray, rhs);
  if(domain->solInfo().order == 2 && C) {
    C->mult(velocity, localTemp);
    rhs.linC(rhs, -1.0, localTemp);
  }
}

void
NonLinDynamic::formRHSpredictor(Vector &velocity, Vector &acceleration, Vector &residual, Vector &rhs, GeomState &geomState, 
                                double midtime, double localDelta)
{
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
      geomState.updatePrescribedDisplacement(userDefineDisp, claw, domain->getNodes());
      setBC(userDefineDisp, userDefineVel);

      // get delta disps
      for(int j = 0; j < claw->numUserDisp; j++)
        userDefineDisp[j] -= userDefineDispLast[j];
 
      // update force residual with KUC
      if(kuc)
        kuc->transposeMultSubtractClaw(userDefineDisp, residual.data(), claw->numUserDisp, clawDofs);

      delete [] userDefineDisp; delete [] userDefineDispLast; delete [] userDefineVel;
    }
  }

  if(domain->solInfo().order == 1) 
    rhs.linC(localDelta, residual);
  else {
    double beta, gamma, alphaf, alpham, dt = 2*localDelta;
    getNewmarkParameters(beta, gamma, alphaf, alpham);
    // rhs = dt*dt*beta*residual + (dt*(1-alpham)*M - dt*dt*(beta-(1-alphaf)*gamma)*C)*velocity
    //       + (dt*dt*((1-alpham)/2-beta)*M - dt*dt*dt*(1-alphaf)*(2*beta-gamma)/2*C)*acceleration
    localTemp.linC(dt*(1-alpham), velocity, dt*dt*((1-alpham)/2-beta), acceleration);
    M->mult(localTemp, rhs);
    if(C) {
      localTemp.linC(-dt*dt*(beta-(1-alphaf)*gamma), velocity, -dt*dt*dt*(1-alphaf)*(2*beta-gamma)/2, acceleration);
      C->multAdd(localTemp, rhs);
    }
    rhs.linAdd(dt*dt*beta, residual);
  }

  times->predictorTime += getTime();
}

double
NonLinDynamic::formRHScorrector(Vector &inc_displacement, Vector &velocity, Vector &acceleration,
                                Vector &residual, Vector &rhs, double localDelta)
{
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
    localTemp.linC(-(1-alpham)/(1-alphaf), inc_displacement, dt*(1-alpham), velocity, dt*dt*((1-alpham)/2-beta), acceleration);
    M->mult(localTemp, rhs);
    if(C) {
      localTemp.linC(-dt*gamma, inc_displacement, -dt*dt*(beta-(1-alphaf)*gamma), velocity, -dt*dt*dt*(1-alphaf)*(2*beta-gamma)/2, acceleration);
      C->multAdd(localTemp, rhs);
    }
    rhs.linAdd(dt*dt*beta, residual);
  }
  times->correctorTime += getTime();
  return rhs.norm();
}

void
NonLinDynamic::processLastOutput()  {

  OutputInfo *oinfo = geoSource->getOutputInfo();
  domain->solInfo().lastIt = true;
  for (int iOut = 0; iOut < geoSource->getNumOutInfo(); iOut++)
    oinfo[iOut].interval = 1;
}

void
NonLinDynamic::buildOps(AllOps<double> & allOps, double Kcoef, double Mcoef, double Ccoef, Rbm * rigidBodyModes)
{
  domain->buildOps<double>(allOps, Kcoef, Mcoef, Ccoef, rigidBodyModes);
}

void
NonLinDynamic::preProcess()
{

 // Allocate space for the Static Timers
 times = new StaticTimers;

 this->openResidualFile();

 totIter = 0;
 fprintf(res,"Iteration Time           Residual\trel. res\tdv\t rel. dv\n");

 // Set the nonlinear tolerance
 tolerance = domain->solInfo().getNLInfo().tolRes;

 // Makes renumbering, connectivities and dofsets
 startTimerMemory(times->preProcess, times->memoryPreProcess);
 domain->preProcessing();

 int numdof = domain->numdof();

 int *bc = (int *) dbg_alloca(sizeof(int)*numdof);
 bcx      = new double[numdof];

 // vcx stores the prescribed velocities, which are associated
 // with prescribed displacements. If a user defined displacement
 // is used, then vcx will contain the user defined velocities also.
 vcx      = new double[numdof];

 int i;
 for(i=0; i<numdof; ++i)
   vcx[i] = 0.0;

 BCond* iVel = domain->getInitVelocity();

 // Make the boundary conditions info
 domain->make_bc( bc, bcx );

 // Now, call make_constrainedDSA(bc) to 
 // built c_dsa that will incorporate all 
 // the boundary conditions info
 domain->make_constrainedDSA(bc);

 // ... SET INITIAL VELOCITY
 for(i = 0; i < domain->numInitVelocity(); ++i) {
   int dof = domain->getCDSA()->locate(iVel[i].nnum, 1 << iVel[i].dofnum);
   if(dof >= 0)
     vcx[dof] = iVel[i].val;
 }

 domain->makeAllDOFs();

 AllOps<double> allOps;

 allOps.M = domain->constructDBSparseMatrix<double>();

 allOps.C = (domain->solInfo().alphaDamp != 0.0 || domain->solInfo().betaDamp != 0.0) ? domain->constructDBSparseMatrix<double>() : 0; // DAMPING

 // TDL
 allOps.Kuc = domain->constructCuCSparse<double>();

 // for initialization step just build M^{-1}
 double Kcoef = 0.0;
 double Mcoef = 1.0;
 double Ccoef = 0.0;

 // HAI
 //int useProjector=domain->solInfo().filterFlags;
 //Rbm *rigidBodyModes = (useProjector || domain->solInfo().rbmflg == 1) ? domain->constructRbm() : 0; // PJSA 9-18-2006

 Rbm *rigidBodyModes = 0;

 int useRbmFilter = domain->solInfo().filterFlags;
 int useGrbm = domain->solInfo().rbmflg;
 int useHzemFilter = domain->solInfo().hzemFilterFlag;
 int useHzem   = domain->solInfo().hzemFlag;

 if (useGrbm || useRbmFilter)
   rigidBodyModes = domain->constructRbm();
 else if(useHzem || useHzemFilter)
   rigidBodyModes = domain->constructHzem();

 buildOps(allOps, Kcoef, Mcoef, Ccoef, (Rbm *) 0); // don't use Rbm's to factor in dynamics

 if(useRbmFilter == 1)
    fprintf(stderr," ... RBM filter Level 1 Requested    ...\n");
 if(useRbmFilter == 2)
    fprintf(stderr," ... RBM filter Level 2 Requested    ...\n");
 if(useHzemFilter)
    fprintf(stderr," ... HZEM filter Requested    ...\n");

 if(useRbmFilter || useHzemFilter)
   projector_prep(rigidBodyModes, allOps.M);

 // TDL Change
 kuc = allOps.Kuc;
 M      = allOps.M;
 C      = allOps.C;
 solver = allOps.sysSolver;
 spm    = allOps.spm;
 prec   = allOps.prec;
 spp    = allOps.spp;

 // ... ALLOCATE MEMORY FOR THE ARRAY OF COROTATORS
 allCorot = new Corotator *[domain->numElements()];

 // ... CREATE THE ARRAY OF POINTERS TO COROTATORS
 domain->createCorotators(allCorot);

 // ... CREATE THE ARRAY OF ELEMENT STIFFNESS MATRICES
 if(C) domain->createKelArray(kelArray, melArray, celArray);
 else domain->createKelArray(kelArray, melArray);

 // Look if there is a user supplied routine for control
 claw = geoSource->getControlLaw();

 // create list of usdd node dofs mapped to cdsa dof numbers
 if(claw)  {
   int nClaw = claw->numUserDisp;
   clawDofs = new int[nClaw];
   for (int j = 0; j < nClaw; ++j) {
     int dd = domain->getDSA()->locate(claw->userDisp[j].nnum, (1 << claw->userDisp[j].dofnum));
     clawDofs[j] = domain->getCDSA()->invRCN(dd);
   }
 }

 // Check to see if there is a user supplied function
 // for displacements, forces or control law
 userSupFunc = domain->getUserSuppliedFunction();
 prevFrc = new PrevFrc(domain->numUncon());

 localTemp.initialize(solVecInfo());

 stopTimerMemory(times->preProcess, times->memoryPreProcess);

}

void NonLinDynamic::openResidualFile()
{
  res = fopen("residuals", "w");
}

Solver *
NonLinDynamic::getSolver()
{
  return solver;
}

SDDynamPostProcessor* 
NonLinDynamic::getPostProcessor()
{
 return new SDDynamPostProcessor( domain, bcx, vcx, times);
}

void
NonLinDynamic::printTimers(double timeLoop)
{

 long memoryUsed = solver->size();
 double solveTime = solver->getSolutionTime();

 times->printStaticTimers( solveTime, memoryUsed, domain, timeLoop );
}

int
NonLinDynamic::getAeroAlg()
{
  return domain->solInfo().aeroFlag;
}

int
NonLinDynamic::getThermoeFlag()
{
  return domain->solInfo().thermoeFlag;
}

int
NonLinDynamic::getThermohFlag()
{
  return domain->solInfo().thermohFlag;
}

int
NonLinDynamic::getAeroheatFlag()
{
  return domain->solInfo().aeroheatFlag;
}


void 
NonLinDynamic::dynamCommToFluid(GeomState* geomState, GeomState* bkGeomState,
                                Vector& velocity, Vector& bkVelocity,
                                Vector& vp, Vector& bkVp, int step, int parity, 
                                int aeroAlg)
{

  times->output -= getTime();

  if(domain->solInfo().aeroFlag >= 0 && !domain->solInfo().lastIt) {
    domain->getTimers().sendFluidTime -= getTime();

    // Make d_n_aero from geomState
    ConstrainedDSA *c_dsa = domain->getCDSA();
    DofSetArray *dsa = domain->getDSA();

    // Note that d_n and a_n are vectors being allocated and de-allocated at
    // each time-step being executed.

    Vector d_n( domain->numUncon(), 0.0 );
    Vector d_n2( domain->numUncon(), 0.0 );

    CoordSet &nodes = domain->getNodes();
    int numNodes = nodes.size();

    for(int i=0; i<numNodes; ++i) {

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

    if(!parity && aeroAlg==5){ 
      for(int i=0; i<domain->numNodes(); ++i) {

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

      d_n.linC(0.5,d_n,0.5,d_n2);
      velocity.linC(0.5,velocity,0.5,bkVelocity);
      vp.linC(0.5,vp,0.5,bkVp);
    }
    Vector a_n( domain->numUncon(), 0.0 );
    

    State state( c_dsa, dsa, bcx, vcx, d_n, velocity, a_n, vp );

    domain->getFileExchanger()->sendDisplacements(state);
    if(verboseFlag) fprintf(stderr," ... Sent displacements to Fluid at step %d\n",(step+1));

    domain->getTimers().sendFluidTime += getTime();
  }

  if(domain->solInfo().aeroheatFlag >= 0) {
    domain->getTimers().sendFluidTime -= getTime();

    // Make d_n_aero from geomState
    ConstrainedDSA *c_dsa = domain->getCDSA();
    DofSetArray *dsa = domain->getDSA();

    // Note that d_n and a_n are vectors being allocated and de-allocated at
    // each time-step being executed.

    Vector d_n( domain->numUncon(), 0.0 );
    Vector d_n2( domain->numUncon(), 0.0 );

    CoordSet &nodes = domain->getNodes();

    for(int i=0; i<domain->numNodes(); ++i) {

      int tloc  = c_dsa->locate(i, DofSet::Temp );
      int tloc1 =   dsa->locate(i, DofSet::Temp );

      if(tloc >= 0)
        d_n[tloc]  = (*geomState)[i].x;
      else if (tloc1 >= 0)
        bcx[tloc1] = (*geomState)[i].x;
    }

    State tempState(c_dsa, dsa, bcx, d_n, velocity, vp);
    domain->getFileExchanger()->sendTemperature(tempState);
    if(verboseFlag) fprintf(stderr," ... [T] Sent temperatures ...\n");

    domain->getTimers().sendFluidTime += getTime();
  }

  if(domain->solInfo().thermohFlag >= 0) {
    /* we have to send the vector of temperatures in NODAL order, not
       in DOF order (in which is d_n)! */

    Vector tempsent(domain->numNodes());

    for(int iNode=0; iNode<domain->numNodes(); ++iNode)
      tempsent[iNode] = (*geomState)[iNode].x;

    domain->getFileExchanger()->sendStrucTemp(tempsent);
    if(verboseFlag) fprintf(stderr," ... [T] Sent temperatures ...\n");
  }

  

  times->output += getTime();

}

void
NonLinDynamic::dynamOutput(GeomState* geomState, Vector& velocity,
                           Vector& vp, double time, int step, Vector& force, 
                           Vector &aeroF, Vector &acceleration) const
{
  times->output -= getTime();

  // PJSA 4-9-08: update geomState for time dependent prescribed displacements and velocities (previously done in computeExternalForce2)
  ControlLawInfo *claw = geoSource->getControlLaw();
  ControlInterface *userSupFunc = domain->getUserSuppliedFunction();
  if(claw && claw->numUserDisp) {
    double *userDefineDisp = new double[claw->numUserDisp];
    double *userDefineVel  = new double[claw->numUserDisp];
    userSupFunc->usd_disp(time,userDefineDisp,userDefineVel);
    DofSetArray *dsa = domain->getDSA();
    for(int i = 0; i < claw->numUserDisp; ++i) {
      int dof = dsa->locate(claw->userDisp[i].nnum,1 << claw->userDisp[i].dofnum);
      if(dof >= 0) {
        bcx[dof] = userDefineDisp[i];  // actually, prescribed displacements are output from the geomState, not bcx
        vcx[dof] = userDefineVel[i];
      }
    }
    geomState->updatePrescribedDisplacement(userDefineDisp, claw, domain->getNodes());
    delete [] userDefineDisp; delete [] userDefineVel;
  }

  domain->postProcessing(geomState, force, aeroF, time, (step+1), velocity.data(), vcx,
                         allCorot, melArray, acceleration.data(), (double *)0 /*acx*/);
  times->output += getTime();
}

void
NonLinDynamic::updatePrescribedDisplacement(GeomState *geomState)
{
 if(domain->solInfo().initialTime == 0.0) {
   // Measure time necessary to update the prescribed displacments
   times->timePresc -= getTime();

   // note 1: "if both IDISP and IDISP6 are present in the input file, FEM selects IDISP6 to initialize the displacement field"
   if((domain->numInitDisp() > 0) && (domain->numInitDisp6() == 0))
     geomState->updatePrescribedDisplacement(domain->getInitDisp(), domain->numInitDisp());
   
   if(domain->numInitDisp6() > 0) 
     geomState->updatePrescribedDisplacement(domain->getInitDisp6(), domain->numInitDisp6());
   
   if(domain->nDirichlet() > 0) 
     geomState->updatePrescribedDisplacement(domain->getDBC(), domain->nDirichlet()); 

   times->timePresc += getTime();
 }
}

void
NonLinDynamic::setBC(double *userDefineDisplacement, double *userDefineVel)
{
  DofSetArray *dsa = domain->getDSA();
  for(int i = 0; i < claw->numUserDisp; ++i) {
    int dof = dsa->locate(claw->userDisp[i].nnum, 1 << claw->userDisp[i].dofnum);
    if(dof >= 0) {
      bcx[dof] = userDefineDisplacement[i];
      vcx[dof] = userDefineVel[i];
    }
  }
}

void
NonLinDynamic::getInitialTime(int &initTimeIndex, double &initTime)
{
 initTimeIndex = domain->solInfo().initialTimeIndex;
 initTime      = domain->solInfo().initialTime;
}

void
NonLinDynamic::getNewmarkParameters(double &beta, double &gamma,
                                    double &alphaf, double &alpham)
{
 beta  = domain->solInfo().newmarkBeta;
 gamma = domain->solInfo().newmarkGamma;
 alphaf = domain->solInfo().newmarkAlphaF;
 alpham = domain->solInfo().newmarkAlphaM;
}

