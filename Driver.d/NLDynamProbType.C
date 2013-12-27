#include <limits>
#include <Timers.d/GetTime.h>
extern int verboseFlag;
extern int totalNewtonIter;
extern int debugFlag;
extern GeoSource * geoSource;

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
  bool failSafe = domain->solInfo().getNLInfo().failsafe;

  if(domain->solInfo().order == 1)
    filePrint(stderr, " ... Implicit Generalized Midpoint Time Integration Scheme: alpha =  %4.2f ...\n", gamma);
  else
    filePrint(stderr, " ... Implicit Newmark Time Integration Scheme: beta = %4.2f, gamma = %4.2f, alphaf = %4.2f, alpham = %4.2f ...\n", beta, gamma, alphaf, alpham);

  // Allocate Vectors to store external force, residual velocity 
  // and mid-point force
  VecType external_force(probDesc->solVecInfo());
  VecType aeroForce(probDesc->solVecInfo()); aeroForce.zero();
  VecType rhs(probDesc->solVecInfo());
  VecType residual(probDesc->solVecInfo());
  VecType totalRes(probDesc->sysVecInfo());
  VecType rhsCopy(probDesc->solVecInfo());

  // Backup variables and states for Aeroelastic using A5 algorithm
  typename StateUpdate::RefState *bkRefState = 0;
  typename StateUpdate::StateIncr *bkStateIncr = 0;
  GeomType *bkGeomState = 0;
  GeomType *bkStepState = 0;
  VecType *bkVelocity_n = 0, *bkAcceleration = 0, *bkV_p = 0, *bkAeroForce = 0;
  if(probDesc->getAeroAlg() == 5 || failSafe) {
    bkVelocity_n = new VecType(probDesc->solVecInfo()); bkVelocity_n->zero();
    bkAcceleration = new VecType(probDesc->solVecInfo()); bkAcceleration->zero();
    bkV_p = new VecType(probDesc->solVecInfo()); bkV_p->zero();
    bkAeroForce = new VecType(probDesc->solVecInfo()); bkAeroForce->zero();
  }
  int parity = 0; // parity for A5 algorithm

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

  probDesc->readRestartFile(displacement, velocity_n, acceleration, v_p, *geomState);
  probDesc->updatePrescribedDisplacement(geomState);
  geomState->setVelocityAndAcceleration(velocity_n,acceleration);

  GeomType *stepState = probDesc->copyGeomState(geomState);
  stateIncr = StateUpdate::initInc(geomState, &residual);
  refState = (domain->solInfo().soltyp == 2) ? 0 : StateUpdate::initRef(geomState);

  if(aeroAlg == 5 || failSafe) {
    bkRefState = StateUpdate::initRef(geomState);
    bkStateIncr = StateUpdate::initInc(geomState, &residual);
    bkGeomState = probDesc->copyGeomState(geomState);
    bkStepState = probDesc->copyGeomState(geomState);
  }

  // Get max number of iterations
  int maxit = probDesc->getMaxit();

  // Get time step size
  double dt0 = probDesc->getDt();

  // Get delta = dt/2
  double delta = probDesc->getDelta();

  double time;
  int step;
  probDesc->getInitialTime(step, time);

  // Evaluate external force at initial time
  // send init. step as -1 so that comm. w/fluid code is avoided
  probDesc->getExternalForce(external_force, constantForce, -1, time, geomState, elementInternalForce, aeroForce, delta);

  // Solve for initial acceleration: a^0 = M^{-1}(fext^0 - fint^0 - C*v^0)
  if(domain->solInfo().iacc_switch && geoSource->getCheckFileInfo()->lastRestartFile == 0) {
    if(domain->solInfo().order == 1) {
      if(verboseFlag) filePrint(stderr," ... Computing initial first time derivative of temperature ...\n");
      probDesc->formRHSinitializer(external_force, velocity_n, elementInternalForce, *geomState, velocity_n);
      probDesc->reBuild(*geomState, 0, delta, time);
      probDesc->getSolver()->reSolve(velocity_n);
      geomState->setVelocity(velocity_n);
    }
    else {
      if(verboseFlag) filePrint(stderr," ... Computing initial acceleration ...\n");
      probDesc->formRHSinitializer(external_force, velocity_n, elementInternalForce, *geomState, acceleration, refState);
      probDesc->reBuild(*geomState, 0, delta, time);
      probDesc->getSolver()->reSolve(acceleration);
      geomState->setAcceleration(acceleration);
    }
  }

  // Output initial geometry state of problem and open output files
  probDesc->dynamOutput(geomState, velocity_n, v_p, time, -1, external_force, aeroForce, acceleration, refState);

  // Get maximum number of iterations
  int maxStep = probDesc->getMaxStep();

  // Begin time marching
  double timeLoop =- getTime();
  double midtime;
  char ch[4] = { '|', '/', '-', '\\' };
  bool useTolInc = (domain->solInfo().getNLInfo().tolInc != std::numeric_limits<double>::infinity() 
                 || domain->solInfo().getNLInfo().absTolInc != std::numeric_limits<double>::infinity());
  double dt = dt0;
  bool failed = false;
  int numConverged = 0;
  int p = 0, q = 1;

  double t0 = time;
  double tmax = maxStep*dt0 + 10*std::numeric_limits<double>::epsilon();

  // Time stepping loop
  for(step = 0; time+dt0/q <= tmax || failed; ) {

    dt = dt0/q;
    domain->solInfo().setTimeStep(dt);
    delta = dt/2;

    if(aeroAlg < 0) {
      filePrint(stderr,"\r  %c  Time Integration Loop: t = %9.3e, dt = %9.3e, %3d%% complete ",
                ch[int((timeLoop + getTime())/250.)%4], time, dt, int((time-t0)/(tmax-t0)*100+0.5));
      if(verboseFlag) filePrint(stderr,"\n");
    }

    if(aeroAlg == 5 || failSafe) {
      if(aeroAlg == 5 && parity == 0 || failSafe && !failed) { // copy current state to backup state
        StateUpdate::copyTo(refState, geomState, stepState, stateIncr, velocity_n, acceleration, v_p, aeroForce,
                            bkRefState, bkGeomState, bkStepState, bkStateIncr, *bkVelocity_n, *bkAcceleration, *bkV_p, *bkAeroForce);
      }
      else { // restore current state to backup state
        StateUpdate::copyTo(bkRefState, bkGeomState, bkStepState, bkStateIncr, *bkVelocity_n, *bkAcceleration, *bkV_p, *bkAeroForce,
                            refState, geomState, stepState, stateIncr, velocity_n, acceleration, v_p, aeroForce);
      }
    }

    double midtime = time + dt*(1 - alphaf);
    
    if(!failed && (aeroAlg < 0 || p%q == 0)) { // for coupled aero, only get fluid load at increments of dt0
      double midtimeExt = (aeroAlg < 0) ? midtime : time + dt0*(1-alphaf);
      double deltaExt = (aeroAlg < 0) ? delta : dt0/2;
      probDesc->getExternalForce(external_force, constantForce, step, midtimeExt, geomState, elementInternalForce, aeroForce, deltaExt);
    }

    double resN, initialRes, err;
    int converged;
    bool feasible;

    // Initialize states
    if(domain->solInfo().soltyp != 2) StateUpdate::copyState(geomState, refState);
    probDesc->initializeParameters(geomState);

    // Constraint enforcement iteration loop
    for(int i = 0; i < domain->solInfo().num_penalty_its; ++i) {

      StateUpdate::zeroInc(stateIncr);

      // Newton iteration loop
      for(int iter = 0; iter < maxit; ++iter) {

        residual = external_force;

        try {         
          // Add stateIncr to geomState and compute element tangent stiffness and internal/follower forces
          StateUpdate::integrate(probDesc, refState, geomState, stateIncr, residual,
                                 elementInternalForce, totalRes, velocity_n,
                                 acceleration, midtime);

          // Compute incremental displacements
          StateUpdate::get_inc_displacement(probDesc, geomState, inc_displac, stepState, domain->solInfo().zeroRot);

          // Form rhs = delta^2*residual - M(inc_displac - delta*velocity_n)
          StateUpdate::formRHScorrector(probDesc, inc_displac, velocity_n,
                                        acceleration, residual, rhs, geomState, delta);

          // in this case of "constraints direct" we can't compute the residual norm until after the solver is rebuilt
          if(domain->solInfo().mpcDirect) probDesc->reBuild(*geomState, iter, delta, midtime);

          // Compute and store the residual norm
          resN = probDesc->getResidualNorm(rhs, *geomState, delta);
          if(iter == 0) initialRes = resN;

          // Copy rhs if necessary before it gets overwritten
          if(probDesc->linesearch().type != 0) rhsCopy = rhs;

          // If the convergence criteria does not involve the solution increment, then 
          // check for convergence now (to avoid potentially unnecessary solve)
          if(useTolInc || !(converged = probDesc->checkConvergence(iter, resN, rhs, rhs, midtime)) ) {

            // Assemble global tangent stiffness
            if(!domain->solInfo().mpcDirect) probDesc->reBuild(*geomState, iter, delta, midtime);

            residual = rhs;
            totalNewtonIter++;

            // Solve ([M] + delta^2 [K])dv = rhs (where rhs is overwritten)
            probDesc->getSolver()->reSolve(rhs);
            probDesc->getConstraintMultipliers(*geomState);

            if(probDesc->linesearch().type != 0) {
              // Optional adjustment of the step length
              StateUpdate::linesearch(probDesc, refState, geomState, stateIncr, rhsCopy, elementInternalForce,
                                      totalRes, midtime, external_force, rhs, velocity_n, acceleration,
                                      inc_displac, delta, domain->solInfo().zeroRot, residual);
            }
            StateUpdate::updateIncr(stateIncr, rhs);  // stateIncr = rhs

            // Update state here if the maximum number of iterations is reached
            if(iter == maxit-1) StateUpdate::updateState(probDesc, geomState, *stateIncr);
          }
          // If the converged criteria does involve the solution increment, then
          // check for convergence now
          if(useTolInc) {
            converged = probDesc->checkConvergence(iter, resN, residual, rhs, midtime);
          }
        }
        catch(std::runtime_error& e) {
          if(!failSafe || debugFlag) cerr << "exception: " << e.what() << endl;
          converged = 0;
          break;
        }

        if(converged == 1)
          break;
        else if(converged == -1 && !failSafe)
          filePrint(stderr," ... Warning, Solution diverging\n");

      } // end of Newton iteration loop

      if(converged == 0 && !failSafe && domain->solInfo().getNLInfo().stepUpdateK != std::numeric_limits<int>::max()) {
        filePrint(stderr,"\r *** WARNING: at time %f Newton solver did not reach convergence after %d iterations"
                         " (residual: initial = %9.3e, final = %9.3e, target = %9.3e)\n", 
                         midtime, maxit, initialRes, resN, probDesc->getTolerance());
      }

      if(failed = (failSafe && converged != 1 && resN > domain->solInfo().getNLInfo().failsafe_tol)) {
        // if a Newton solve fails to converge, terminate constraint enforcement iterations
        break;
      }

      // update lagrange multipliers and/or penalty parameters 
      probDesc->updateParameters(geomState);

      // check constraint violation error
      feasible = probDesc->checkConstraintViolation(err);
      if(feasible) break;

    } // end of constraint enforcement iteration loop

    if(failed) { 
      // if a Newton solve fails to converge, decrease the time step and try again
      p *= 2;
      q *= 2;
      numConverged = 0;
      continue;
    }

    // Step Update: updates position from _{n+1-alphaf} to _{n+1} and velocity/acceleration from _{n} to _{n+1}
    v_p = velocity_n;
    StateUpdate::midpointIntegrate(probDesc, velocity_n, delta,
                                   stepState, geomState, stateIncr, residual,
                                   elementInternalForce, totalRes, acceleration,
                                   domain->solInfo().zeroRot);
    if(domain->solInfo().soltyp != 2) {
      probDesc->updateStates(refState, *geomState); // update internal states to _{n+1}
      StateUpdate::copyState(geomState, refState);
    }

    p++;
    numConverged++;
    if(numConverged >= 2 && (q-p)%2 == 0 && q >= 2) { p /= 2; q /= 2; numConverged = 0; } // increase the timestep
    time += dt;

    if(p%q == 0) { // finished a whole step. post process

      // Output results at current time
      if(step+1 == maxStep && (aeroAlg != 5 || parity == 1)) probDesc->processLastOutput();
      else if(aeroAlg >= 0 || probDesc->getThermohFlag() >= 0 || probDesc->getAeroheatFlag() >= 0) {
        probDesc->dynamCommToFluid(geomState, bkGeomState, velocity_n, *bkVelocity_n, v_p, *bkV_p, step, parity, aeroAlg);
      }
      probDesc->dynamOutput(geomState, velocity_n, v_p, time, step, external_force, aeroForce, acceleration, refState);

      if(aeroAlg == 5) { 
        if(!parity) {
          time -= dt;
          step--;
        }
        parity = ( parity ? 0 : 1 );
      }

      step = p/q;
    }
  } // end of time stepping loop

  if(aeroAlg < 0)
    filePrint(stderr,"\r ... Time Integration Loop: t = %9.3e, dt = %9.3e, 100%% complete ...\n", time, dt);

  timeLoop += getTime();
#ifdef PRINT_TIMERS
  filePrint(stderr," ... Total Loop Time = %.2e s   ...\n",timeLoop/1000.0);
#endif

  probDesc->printTimers(timeLoop);

  delete stepState;
  delete geomState;

  if(aeroAlg == 5 || failSafe) {
    delete bkVelocity_n;
    delete bkAcceleration;
    delete bkV_p;
    delete bkAeroForce;
    delete bkRefState;
    delete bkStateIncr;
    delete bkGeomState;
    delete bkStepState;
  }

}

