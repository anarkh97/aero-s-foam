#include <Timers.d/GetTime.h>
extern int verboseFlag;
extern int totalNewtonIter;

/****************************************************************
 *
 *  Purpose: Implicit Nonlinear Dynamic Algorithm based on 
 *           generalized-alpha method (template version, to be used in 
 *           single and multiple domain FEM problems) 
 *
 *  Input: Solver type, Vector Type, Post Processor type, 
 *         Problem descriptor type and Geometry State type template arguments
 *
 *  Output:
 *         Geometric state of structure
 *
 *  Coded by: Kendall Pierson and Teymour Manzouri
 *
 *  Date: March 1998
 *
 * Single Domain Template arguments:
 *   OpSolver          = Solver 
 *   VecType           = Vector
 *   PostProcessor     = SDPostProcessor
 *   ProblemDescriptor = NonLinDynamic
 *   GeomType          = GeomState
 ***************************************************************/

#include <Timers.d/GetTime.h>
#include <Corotational.d/GeomState.h>

template < class OpSolver, 
           class VecType, 
           class PostProcessor, 
           class ProblemDescriptor, 
           class GeomType,
           class StateUpdate >
void
NLDynamSolver < OpSolver, VecType, PostProcessor, ProblemDescriptor,
		GeomType, StateUpdate >::solve() 
{
  // Set up nonlinear dynamics problem descriptor 
  probDesc->preProcess();

  // Compute time integration values: dt, totalTime, maxStep
  probDesc->computeTimeInfo();
  probDesc->getNewmarkParameters(beta, gamma, alphaf, alpham);

  // Get pointer to Solver
  OpSolver *solver = probDesc->getSolver();

  if(domain->solInfo().order == 1)
    filePrint(stderr, " ... Implicit Newmark Algorithm     ...\n");
  else
    filePrint(stderr, " ... Implicit Newmark Time Integration Scheme: beta = %4.2f, gamma = %4.2f, alphaf = %4.2f, alpham = %4.2f ...\n", beta, gamma, alphaf, alpham);

  // Allocate Vectors to store external force, residual velocity 
  // and mid-point force
  VecType external_force(probDesc->solVecInfo());
  VecType aeroForce(probDesc->solVecInfo());
  VecType rhs(probDesc->solVecInfo());
  VecType residual(probDesc->solVecInfo());
  VecType totalRes(probDesc->sysVecInfo());

  // Backup variables and states for Aeroelastic using A5 algorithm
  typename StateUpdate::RefState *bkRefState = 0;
  typename StateUpdate::StateIncr *bkStateIncr = 0;
  GeomType *bkGeomState = 0;
  GeomType *bkStepState = 0;
  VecType *bkVelocity_n = 0, *bkAcceleration = 0, *bkV_p = 0, *bkAeroForce = 0;
  if(probDesc->getAeroAlg() == 5) {
    bkVelocity_n = new VecType(probDesc->solVecInfo()); bkVelocity_n->zero();
    bkAcceleration = new VecType(probDesc->solVecInfo()); bkAcceleration->zero();
    bkV_p = new VecType(probDesc->solVecInfo()); bkV_p->zero();
    bkAeroForce = new VecType(probDesc->solVecInfo()); bkAeroForce->zero();
  }
  int parity = 0; //parity for A5 algorithm

  // zero vectors
  external_force.zero();
  rhs.zero();
  residual.zero();
  totalRes.zero();

  VecType elementInternalForce(probDesc->elemVecInfo());
  elementInternalForce.zero();

  VecType inc_displac(probDesc->solVecInfo());
  inc_displac.zero();

  VecType constantForce(probDesc->solVecInfo());
  probDesc->getConstForce(constantForce);

  VecType displacement(probDesc->solVecInfo()); displacement.zero(); 
  VecType velocity_n(probDesc->solVecInfo()); velocity_n.zero();   // velocity at time t_n
  VecType acceleration(probDesc->solVecInfo()); acceleration.zero();
  VecType v_p(probDesc->solVecInfo()); v_p.zero();

  // Set up initial conditions and check for AEROELASTIC computation
  int aeroAlg = probDesc->getInitState(displacement,velocity_n,acceleration,v_p);
  if(aeroAlg == 1 || aeroAlg == 8) return;

  // Initialize geometric state of problem using the mesh geometry,
  // restart file (if it exists), or the initial displacements (if any).
  GeomType *geomState = probDesc->createGeomState();

  stateIncr = StateUpdate::initInc(geomState, &residual);
  refState = StateUpdate::initRef(geomState);

  if(aeroAlg == 5) {
    bkRefState = StateUpdate::initRef(geomState);
    bkStateIncr = StateUpdate::initInc(geomState, &residual);
    bkGeomState = probDesc->copyGeomState(geomState);
    bkStepState = probDesc->copyGeomState(geomState);
  }

  probDesc->readRestartFile(displacement, velocity_n, acceleration, v_p, *geomState);
  probDesc->updatePrescribedDisplacement(geomState);
  GeomType *stepState = probDesc->copyGeomState(geomState);

  // Get max number of iterations
  int maxit = probDesc->getMaxit();

  // Get time step size
  double dt = probDesc->getDt();

  // Get delta = dt/2
  double delta = probDesc->getDelta();

  double time;
  int step;
  probDesc->getInitialTime(step, time);

  // Evaluate external force at initial time
  // send init. step as -1 so that comm. w/fluid code is avoided
  probDesc->getExternalForce(external_force, constantForce, -1, time, geomState, elementInternalForce, aeroForce);

  // Solve for initial acceleration: a^0 = M^{-1}(fext^0 - fint^0 - C*v^0)
  if(domain->solInfo().iacc_switch) {
    if(domain->solInfo().order == 1) {
      if(verboseFlag) filePrint(stderr," ... Computing initial first time derivative of temperature ...\n");
      probDesc->formRHSinitializer(external_force, velocity_n, elementInternalForce, *geomState, velocity_n);
      solver->reSolve(velocity_n);
    }
    else {
      if(verboseFlag) filePrint(stderr," ... Computing initial acceleration ...\n");
      probDesc->formRHSinitializer(external_force, velocity_n, elementInternalForce, *geomState, acceleration);
      solver->reSolve(acceleration);
    }
  }

  // Output initial geometry state of problem and open output files
  probDesc->dynamOutput(geomState, velocity_n, v_p, time, -1, external_force, aeroForce, acceleration);

  // Get maximum number of iterations
  int maxStep = probDesc->getMaxStep();

  // Begin time marching
  double timeLoop =- getTime();
  double currentRes;
  char ch[4] = { '|', '/', '-', '\\' };

  for( ; step < maxStep; ++step) {

    filePrint(stderr,"\r  %c  Time Integration Loop: t = %9.3e, %3d%% complete ",ch[int((timeLoop + getTime())/250.)%4], time+dt, int(double(step+1)/double(maxStep)*100.0));
    if(verboseFlag) cerr << endl;

    if(aeroAlg == 5) {
      if(parity==0) //copy current state to backup state
        StateUpdate::copyTo(refState, geomState, stepState, stateIncr, velocity_n, acceleration, v_p, aeroForce,
                            bkRefState, bkGeomState, bkStepState, bkStateIncr, *bkVelocity_n, *bkAcceleration, *bkV_p, *bkAeroForce);
      else //restore backup state to current state
        StateUpdate::copyTo(bkRefState, bkGeomState, bkStepState, bkStateIncr, *bkVelocity_n, *bkAcceleration, *bkV_p, *bkAeroForce,
                            refState, geomState, stepState, stateIncr, velocity_n, acceleration, v_p, aeroForce);
    }

    double midtime = time + dt*(1 - alphaf);
    probDesc->getExternalForce(external_force, constantForce, step, midtime, geomState, elementInternalForce, aeroForce);

    time += dt;

    double resN, initialRes;
    int converged;

    // Initialize states
    StateUpdate::copyState(geomState, refState);
    StateUpdate::zeroInc(stateIncr);

    // Iteration loop
    for(int iter = 0; iter < maxit; ++iter, ++totalNewtonIter) {

      residual = external_force;
         
      // And stateIncr to geomState and compute element tangent stiffness and internal/follower forces
      // also update the constraints (updateContactConditions called)
      StateUpdate::integrate(probDesc, refState, geomState, stateIncr, residual,
                             elementInternalForce, totalRes, velocity_n,
                             acceleration, midtime);

      // Assemble global tangent stiffness
      probDesc->reBuild(*geomState, iter);

      // Compute incremental displacements
      geomState->get_inc_displacement(inc_displac, *stepState);

      // Form rhs = delta^2*residual - M(inc_displac - delta*velocity_n)
      resN = StateUpdate::formRHScorrector(probDesc, inc_displac, velocity_n,
                                           acceleration, residual, rhs, geomState);
      resN = probDesc->getResidualNorm(rhs); // addMpcForces called

      //filePrint(stderr,"2 NORMS: fext*fext %e %e residual*residual %e\n", external_force*external_force, resN*resN);

      currentRes = resN;
      if(iter == 0) initialRes = resN;
      residual = rhs;

      // Solve ([M] + delta^2 [K])dv = rhs (where rhs is over written)
      //cerr << "rhs*rhs = " << rhs*rhs << endl;
      solver->reSolve(rhs);
      //cerr << "sol*sol = " << rhs*rhs << endl;

      // Check for convergence
      // XXXX it seems like a waste of one rebuild/solve to compute dv before checking for convergence. dv is only used for printing
      converged = probDesc->checkConvergence(iter, resN, residual, rhs, time);
      StateUpdate::updateIncr(stateIncr, rhs);  // stateIncr = rhs

      if(converged == 1)
        break;
      else if(converged == -1)
        filePrint(stderr," ... Warning, Solution diverging\n");
    }
    if(converged == 0) 
      filePrint(stderr,"\r *** WARNING: at time %f Newton solver did not reach convergence after %d iterations (residual: initial = %9.3e, final = %9.3e, target = %9.3e)\n", 
                time, maxit, initialRes, currentRes, probDesc->getTolerance());

    StateUpdate::copyState(geomState, refState);

    // Step Update (updates state which includes displacement and velocity, but not acceleration) 
    v_p = velocity_n;
    StateUpdate::midpointIntegrate(probDesc, velocity_n, delta,
                                   stepState, geomState, stateIncr, residual,
                                   elementInternalForce, totalRes, acceleration); // note: stateIncr is not used in this function except for the TotalUpdater

    // Update the acceleration: a^{n+1} = (v^{n+1}-v^n)/delta - a^n
    if(domain->solInfo().order != 1)
      acceleration.linC(-(1-gamma)/gamma, acceleration, -1/(2*delta*gamma), v_p, 1/(2*delta*gamma), velocity_n);

    // Output results at current time
    if(step+1 == maxStep && (aeroAlg != 5 || parity==1)) probDesc->processLastOutput();
    if(aeroAlg >= 0 || probDesc->getThermohFlag() >= 0 || probDesc->getAeroheatFlag() >= 0) {
      probDesc->dynamCommToFluid(geomState, bkGeomState, velocity_n, *bkVelocity_n, v_p, *bkV_p, step, parity, aeroAlg);
    }
    probDesc->dynamOutput(geomState, velocity_n, v_p, time, step, external_force, aeroForce, acceleration);

    if(aeroAlg == 5) { 
      if(!parity) {
        time -= dt;
        step--;
      }
      parity = ( parity ? 0 : 1 );
    }
  }
  if(!aeroAlg)
    filePrint(stderr,"\r ... Time Integration Loop: t = %9.3e, 100%% complete ...\n", time);

  timeLoop += getTime();
#ifdef PRINT_TIMERS
  filePrint(stderr," ... Total Loop Time = %.2e s   ...\n",timeLoop/1000.0);
#endif

  probDesc->printTimers(timeLoop);
}

