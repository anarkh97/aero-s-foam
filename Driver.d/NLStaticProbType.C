#include <stdio.h>
#include <Timers.d/GetTime.h>

#define DEBUG_NEWTON // use this to output at every newton iteration

extern int verboseFlag;

template < class OpSolver, 
           class VecType, 
	   class PostProcessor, 
           class ProblemDescriptor, 
           class GeomType,
	   class StateUpdate >
void
NLStaticSolver < OpSolver, VecType, PostProcessor, ProblemDescriptor, GeomType, StateUpdate >::solve()
{
 // Set up nonlinear Problem 
 probDesc->preProcess();

 // Get Solver
 OpSolver *solver = probDesc->getSolver();

 // Allocate Appropriate Vectors to store external force and residual
 VecType force(probDesc->solVecInfo());
 VecType residual(probDesc->solVecInfo());
 VecType totalRes(probDesc->sysVecInfo());

 VecType elementInternalForce(probDesc->elemVecInfo());

 residual.zero();
 totalRes.zero();
 elementInternalForce.zero();

 // Get right hand side (external force)
 probDesc->getRHS(force);

 // Initialize geometric state of Problem
 GeomType *geomState = probDesc->createGeomState();
 stateIncr = StateUpdate::initInc(geomState, &residual);
 
 refState = StateUpdate::initRef(geomState);

 double lambda = 0.0;

 // Incremental step parameter for load / prescribed displacement control 
 double deltaLambda = probDesc->getDeltaLambda0(); 

 // Newton iteration loop
 double maxLambda = probDesc->getMaxLambda();

 // Output structure initial configuration
 if(deltaLambda != maxLambda)
   probDesc->staticOutput(geomState, lambda, force, totalRes);

 int numIter = 0;

 int numSteps = int(0.5+(maxLambda / deltaLambda));

 // Set initial value of lambda equal to deltaLambda
 lambda = deltaLambda;

 int step;
 for(step = 1; step <= numSteps; ++step) {

   filePrint(stderr," --------------------------------------\n");
   filePrint(stderr," ... Newton : Start Step #%d --- Lambda = %e\n",step, lambda);

   StateUpdate::copyState(geomState, refState);
   probDesc->updatePrescribedDisplacement(geomState, lambda);


   // call newton iteration with load step lambda
   int converged = newton(force, residual, totalRes,
                          elementInternalForce, solver, 
                          refState, geomState, numIter, lambda, step);

   double time = lambda;
   if(deltaLambda == maxLambda) time = 0.0;

   if(converged == 1)
     filePrint(stderr," ... Newton : Step #%d converged after %d iterations\n",
                    step,numIter+1);
   else if(converged == -1) {
     filePrint(stderr," ... Newton : Step #%d diverged after %d iterations\n",step,
                   numIter);
     filePrint(stderr," ... Newton : analysis interrupted by divergence\n");
#ifndef DEBUG_NEWTON
     probDesc->staticOutput( geomState, time, force, totalRes);
#endif
     break;
   } 

   filePrint(stderr, " ... Newton : End   Step #%d --- Max Steps = %d\n", step, numSteps);

   if(converged ==0)
     filePrint(stderr," *** WARNING: Newton step did not converge after %d iterations (res = %e, target = %e)\n",
               numIter, totalRes.norm(), probDesc->getTolerance());

   fflush(stderr);

   filePrint(stderr," --------------------------------------\n");
   fflush(stderr);
#ifndef DEBUG_NEWTON
   // Output current load step results
   probDesc->staticOutput(geomState, time, force, totalRes);
#endif

   // increment load parameter
   lambda += deltaLambda;
 }

 probDesc->printTimers();

}


template < class OpSolver,
           class VecType,
           class PostProcessor,
           class ProblemDescriptor,
           class GeomType,
	   class StateUpdate>
void
NLStaticSolver < OpSolver, VecType, PostProcessor, ProblemDescriptor, GeomType, StateUpdate >
::arclength()
{
 // WARNING THE CREATION OF REFSTATE MUST BE DONE IF THIS IS TO BE USED
 // WITH TOTALUPDATE MODE
 // Set up  nonlinear Problem
 probDesc->preProcess();

 // Get Solver
 OpSolver *solver = probDesc->getSolver();

 // Allocate Vectors
 VecType force(probDesc->solVecInfo());
 VecType residual(probDesc->solVecInfo());
 VecType totRes(probDesc->sysVecInfo());
 VecType arcLenResid(probDesc->solVecInfo());
 VecType pVec(probDesc->solVecInfo());

 VecType elementInternalForce(probDesc->elemVecInfo());

 // HB
 residual.zero();
 totRes.zero();
 elementInternalForce.zero();

 // Get right hand side
 probDesc->getRHS(force);

 // Compute initial force norm
 double forceNorm = force.norm();

 // Initialize geometric state of Problem
 GeomType *u0 = probDesc->createGeomState();
 GeomType *u  = probDesc->createGeomState();

 // HB
 GeomType *geomState = probDesc->createGeomState();
 stateIncr = StateUpdate::initInc(geomState, &residual);

 refState = StateUpdate::initRef(geomState);

 // ... Output initial configuration
 probDesc->staticOutput(u0, 0.0, force, totRes);

 // ... DEFINE deltaLambda0
 double deltaLambda0 = probDesc->getDeltaLambda0();

 // ... DEFINE minimum delta Lambda and maximum delta Lambda
 // lmin <= deltaLamba <= lmax
 double lfactor = probDesc->getScaleFactor();
 double lmin = deltaLambda0 / lfactor;
 double lmax = deltaLambda0 * lfactor;

 // ... DEFINE initial lambda0
 double lambda0 = 0.0;

 // ... COMPUTE FIRST LOAD CONTROL PARAMETER
 double lambda = lambda0 + deltaLambda0;

 int numIter = 0;
 int step    = 1;
 newton(force, residual, elementInternalForce, totRes, solver, refState,
        u, numIter, lambda, step);

 // ... Declare Vector dU
 VecType dU(probDesc->solVecInfo());

 // ... Define initial dU = u - u0
 u->diff(*u0, dU);

 // ... Define deltaS   (arc length distance)
 double deltaS = dU.norm();
 filePrint(stderr," -> deltaS = %e\n",deltaS);
 double deltaSmax = 10*deltaS;
 double deltaSmin = 0.1*deltaS;

 // ... COMPUTE W = SMOOTHING PARAMETER
 //double w = (deltaS*deltaS) / (deltaLambda0*deltaLambda0);
 //HB: the same as in Corotational.d/GeomNLSolver.C arcLength()
 // -> it seems to work but for me this is wrong !!!
 //deltaS = deltaS*deltaS;
 double w = (deltaS*deltaS) / (deltaLambda0*deltaLambda0);
 w = 0;

 // ... DEFINE NU = load control parameter multiplier
 double nu = 1.0;

 // ... DEFINE MAXIMUM NUMBER OF TRAJECTORY ITERATIONS
 int maxNumTrajectory = 2000;
 //HB
 //maxNumTrajectory = 1;

// ---------------------------------------------------

 double deltaLambda;
 int iter, numExtIter;

  // ... COMPUTE TRAJECTORY (EQUILIBRIUM) PATH
 for(iter=0; iter<maxNumTrajectory; ++iter) {

   probDesc->staticOutput(u, lambda, force, totRes);

   // Compute: dU = nu*(u - u0);
   u->diff(*u0, dU);
   filePrint(stderr," -> ||dU|| = %e\n",dU.norm());
   //dU *= nu;

   *u0 = *u;

   // ... COMPUTE NEXT DELTALAMBDA
   //deltaLambda = nu*(lambda - lambda0);
   deltaLambda = (lambda - lambda0);

   // ... CHECK MAGNITUDE OF DELTALAMBDA
   if(deltaLambda > lmax ) deltaLambda = lmax;
   if(deltaLambda < lmin ) deltaLambda = lmin;

   // ... SET lambda0 = lambda
   lambda0 = lambda;

   // ... COMPUTE DELTAS
   double deltaS0 = deltaS;
   deltaS *= nu;
   if(deltaS > deltaSmax) deltaS = deltaSmax;
   if(deltaS < deltaSmin) deltaS = deltaSmin;

   filePrint(stderr,"**** Equilibrium #%d deltaLambda = %10.6f lambda = %10.6f\n"
                   ,step, deltaLambda, lambda);
   step++;

   predictorStep(*u, *u0, dU, lambda, deltaLambda, deltaS, deltaS0, w, solver,
                 force, residual,totRes, elementInternalForce,  pVec, step);

   //u->diff(*u0, dU);
   filePrint(stderr," -> ||dU|| = %e\n",dU.norm());

   // ... CALL EXTENDED NEWTON
   extendedNewton(*u,dU,lambda,deltaLambda,deltaS,w,numExtIter,
                  solver,force,residual, totRes, arcLenResid, forceNorm,
                  elementInternalForce, pVec, step);

   if( abs(lambda) >= 1.0 ) break;

   //...DETERMINE CONTROL PARAMETER BASED ON # OF ITERATIONS IN EXTENDED NEWTON
   nu = 1.0;
   //HB: disable this for debug
   //nu = sqrt(4./numExtIter);
   nu = pow(4./numExtIter,0.75);
   //if(numExtIter < 4) nu = 2.0;
   //if(numExtIter > 6) nu = 0.5;

 }
 //HB
 filePrint(stderr,"**** Equilibrium #%d deltaLambda = %10.6f lambda = %10.6f\n"
                 ,step, deltaLambda, lambda);


 // ... DEFINE INITIAL GUESS
 u->interp((1.0 - lambda0)/(lambda - lambda0), *u, *u0 );
 lambda = 1.0;

 // ... CALL NEWTON FOR FINAL SOLUTION
 step++;
 newton(force, residual, totRes, elementInternalForce, solver, refState, u, numIter, lambda, step);

 // CALL POST PROCESSING OF DISPLACEMENTS
 probDesc->staticOutput(u, lambda, force, totRes);
 
}


template < class OpSolver,
           class VecType,
           class PostProcessor,
           class ProblemDescriptor,
           class GeomType,
	   class StateUpdate>
int
NLStaticSolver < OpSolver, VecType, PostProcessor, ProblemDescriptor, GeomType, StateUpdate  >
::newton( VecType& force, VecType& residual, VecType &totalRes,
          VecType& elementInternalForce, 
          OpSolver* solver, typename StateUpdate::RefState *refState,
	  GeomType* geomState,
	  int &numIter, double lambda, 
          int step )
{
  // Accumulate time spent in solving and geomstate update for one step
  double timeSolve   = 0.0;
  double timeStiff   = 0.0;
  double timeRebuild = 0.0;

  int maxit = probDesc->getMaxit();
  
  // Zero the state increment
  StateUpdate::zeroInc(stateIncr);
  
  // Main Newton Iteration Loop
  double e_k;
  int iter, converged;
  for(iter = 0; iter < maxit; ++iter) {

    // residual = lambda*force;
    residual.linC(force, lambda);
 
    // Update geomState then compute current tangent stiffness and residual force (including follower force contributions)
    timeStiff -= getTime();
    double residualNorm = StateUpdate::integrate(probDesc, refState, geomState, stateIncr,
                                                 residual, elementInternalForce, totalRes, lambda);
    e_k = probDesc->getEnergy(lambda, force, geomState);
    timeStiff += getTime();
#ifdef PRINT_TIMERS
    filePrint(stderr,"  Rebuild Element Stiffness & Internal Force time = %13.4f s\n", 
              timeStiff/1000.0);
#endif    

    // Rebuild tangent stiffness matrix when necessary
    timeRebuild -= getTime();
    int rebuildFlag = probDesc->reBuild(iter, step, *geomState);
    timeRebuild += getTime();
#ifdef PRINT_TIMERS
    filePrint(stderr,"  Rebuild Tangent Stiffness Matrix time = %13.4f s\n",
              timeRebuild/1000.0);
#endif

    if(rebuildFlag) {
      filePrint(stderr," ... Newton : Iter #%d --- Rebuild Tangent Stiffness (res = %e)\n", iter+1, residualNorm); // HB
    }
    else filePrint(stderr," ... Newton : Iter #%d (res = %e)\n", iter+1, residualNorm);

    // Solve current system Kt*u = residual, overwrite residual with u
    timeSolve -= getTime();
    solver->reSolve(residual);
    timeSolve += getTime();
#ifdef PRINT_TIMERS
    filePrint(stderr,"  Solve Incremental Displacement %13.4f s\n",
              timeSolve/1000.0);
#endif

    if(probDesc->linesearch()) { // experimental
      double alpha, alpha_opt = std::numeric_limits<double>::epsilon();
      VecType tmp(probDesc->solVecInfo());
      for(alpha = 1.0e3; alpha > 0.001; alpha *= 0.9) {
        GeomType *tmpState = new GeomType(*geomState);
        tmp.linC(residual,alpha);
        StateUpdate::updateIncr(stateIncr, tmp);
        StateUpdate::integrate(probDesc, refState, tmpState, stateIncr, tmp, elementInternalForce, totalRes);
        double e = probDesc->getEnergy(lambda, force, tmpState);
        cerr << "alpha = " << alpha << ", e = " << e << endl;
        delete tmpState;
        if(e < e_k) { e_k = e; alpha_opt = alpha; }
      }
      cerr << "alpha_opt = " << alpha_opt << endl;
      residual *= alpha_opt;
    }

    StateUpdate::updateIncr(stateIncr, residual);

    // Compute incremental displacement norm
    double normDv = residual.norm();

    // Check convergence using residual norm & incremental displacement norm
    converged = probDesc->checkConvergence(iter, normDv, residualNorm);

#ifdef DEBUG_NEWTON
    probDesc->staticOutput( geomState, double(iter), force, totalRes);
#endif

    // If converged, break out of loop
    if(converged == 1) break; // PJSA_DEBUG don't test for divergence
  }

  // return with the number of iterations newton took to converge/diverge
  numIter = iter;

  return converged;
}

// HB
template < class OpSolver,
           class VecType,
           class PostProcessor,
           class ProblemDescriptor,
           class GeomType,
           class StateUpdate>
void
NLStaticSolver < OpSolver, VecType, PostProcessor, ProblemDescriptor, GeomType, StateUpdate >
::predictorStep(GeomType &u, GeomType &u0, VecType &dU, double &lambda, double &deltaLambda,
                 double &deltaS, double &deltaS0, double w, OpSolver* solver,
                 VecType& force, VecType& residual, VecType &totRes, VecType& elementInternalForce, 
                 VecType& duds, int step)
{
  filePrint(stderr," ### GET IN PREDICTOR STEP\n");
  //HB: be aware that for multi-domain problem most of the dot product here are not "consistent"
  //    in the sense that they won't give the same answer as for the same non multi-domains problem
  //    (see the way a dot product is computed for a GenDistrVector). But in the case were the 
  //    distributed vector has the same value on its share nodes, the computed dot product can be
  //    interpreted as a scaled/weighted dot product: this is why it should be ok to use it here
  //    But be carefull, if we use this with a residual-like vector (i.e. a "disassembled-like vector") ...
  probDesc->getRHS(force); //PJSA nonlinear external force computed in getStiffAndForce probDesc->getRHS(force, &u);
  residual.zero();
  double residualNorm = probDesc->getStiffAndForce(u, residual, elementInternalForce, totRes, lambda);
  int rebuildFlag = probDesc->reBuild(0, step, u);
  if(rebuildFlag) { filePrint(stderr," ... ||res|| = %e\n", residualNorm); }

  duds.zero(); 
  solver->solve(force, duds);

  double m = 2*(dU*duds);
  double n = 2*w*deltaLambda*force.sqNorm();
  //double dlds = 1./(m+n);
  double dlds = 2.*deltaS0/(m+n);
  filePrint(stderr," -> m = %e, n = %e, dlds = %e\n",m,n,dlds);
  duds *= (dlds*deltaS);
  filePrint(stderr," -> ||duds||    = %e\n",duds.norm());

  u.update(duds);
  deltaLambda = dlds*deltaS;
  lambda += deltaLambda;
  filePrint(stderr," -> deltaLambda = %e\n",deltaLambda);
  u.diff(u0,dU);
  filePrint(stderr," -> ||dU||      = %e\n",dU.norm());
  //dU = duds;
}

template < class OpSolver,
           class VecType,
           class PostProcessor,
           class ProblemDescriptor,
           class GeomType,
	   class StateUpdate>
void
NLStaticSolver < OpSolver, VecType, PostProcessor, ProblemDescriptor, GeomType, StateUpdate >
::extendedNewton(GeomType &u, VecType &dU, double &lambda, double deltaLambda, 
                 double &deltaS, double w, int &numExtIter, OpSolver* solver,
                 VecType& force, VecType& residual,  VecType &totRes,
                 VecType& arcLenResid, 
                 double forceNorm, VecType& elementInternalForce, VecType& pVec, int step)
{

  //HB
  filePrint(stderr," #####################################\n");
  filePrint(stderr," ### GET IN EXTENDED-NEWTON SOLVER ###\n");
  filePrint(stderr," #####################################\n");
  // Compute new lambda
  double lambda0 = lambda-deltaLambda;
  //lambda += deltaLambda;

  // Compute first norm
  double firstNorm = (lambda*lambda)*forceNorm;
  //HB
  fprintf(stderr,"--- firstNorm = %e\n",firstNorm);


  if(firstNorm == 0.0) return;

  //predictorStep(u, lambda, deltaLambda, deltaS, w, solver,
  //              force, residual, pVec, step);


  //double Dlbda = lambda-lambda0;
  double Dlbda = deltaLambda;
  double r0= (dU*dU) + w*(Dlbda*Dlbda) - (deltaS*deltaS);
  //r0 = 1.;

  int maxExtIter = probDesc->getMaxit();
  for(numExtIter = 0; numExtIter < maxExtIter; ++numExtIter) {

    filePrint(stderr," ------------------------------------\n");
    filePrint(stderr," ### Extended-Newton iteration %d ###\n",numExtIter);
    filePrint(stderr," ------------------------------------\n");
    // HB 
    //probDesc->getRHS(force, &u);

    // Compute residual = lambda*force
    residual.linC(force, lambda);

    // Compute stiffness and residual force
    double residualNorm = probDesc->getStiffAndForce(u, residual, elementInternalForce, totRes, lambda);

    // rebuild tangent stiffness matrix when necessary
    // step # should be the second argument in probDesc->reBuild() 
    //probDesc->reBuild(numExtIter, 1, u);
    int rebuildFlag = probDesc->reBuild(numExtIter, step, u); 

    if(rebuildFlag) { filePrint(stderr," ... ||res|| = %e\n", residualNorm); }

    //HB
    Dlbda = lambda-lambda0;
    double r = (dU*dU) + w*(Dlbda*Dlbda) - (deltaS*deltaS);
    fprintf(stderr," r = %e, |r/r0| = %e,|dl| = %e, ||dU|| = %e\n",r, fabs(r/r0), fabs(Dlbda), sqrt(dU*dU));

    arcLenResid = residual;

    filePrint(stderr," ### First solve ###\n");
    solver->reSolve(residual);

    //solver->resetOrthoSet(); // HB: force reseting orthoset 
    filePrint(stderr," ### Second solve ###\n");
    solver->solve(force, pVec);

    //HB
    //double mu = -1.0*( dU * residual ) / (w*deltaLambda + dU * pVec);
    double mu = -0.5*( r + 2.*(dU*residual)) / (w*Dlbda + dU*pVec);

    // Compute: pVec = mu*pVec + residual
    pVec.linC(mu, pVec, 1.0, residual);
    //HB
    dU.linC(dU,1.0,pVec);
    //u.diff(*u0,dU);

    // Update geometrical state
    u.update(pVec);

    // Update lambda
    lambda += mu;
    // HB
    fprintf(stderr,"              lambda = %f, mu = %15.6e\n",lambda,mu);

    // Compute arcLenResid = arcLenResid + mu*force
    arcLenResid.linC(arcLenResid, mu, force);

    // double residualNorm = arcLenResid.norm();
    double normDv       = pVec.norm();

    // KHP: same problem as before, norm of force contains 
    //      lagrange multiplier forces and thus norm does not go to zero.
    //int converged = probDesc->checkConvergence(numExtIter, normDv, 
    //                residual.norm());
    // HB
    int converged = probDesc->checkConvergence(numExtIter, normDv, residualNorm);

    if (converged) break;
  }

  filePrint(stderr,"... Number of External Iterations = %d\n", numExtIter);

}
