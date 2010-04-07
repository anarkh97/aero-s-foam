#include <Pita.d/NLDynamTimeIntegrator.h>
#include <Problems.d/NonLinDynam.h>
#include <Solvers.d/Solver.h>
#include <cstdio>
#include <cstdlib>

using std::fprintf;
using std::exit;

extern int verboseFlag;

// Constructor & destructor
// ------------------------
NLDynamTimeIntegrator::NLDynamTimeIntegrator(NonLinDynamic & pbDesc) :
  probDesc(pbDesc),
  postProcessor_(&(pbDesc.defaultPostProcessor())),
  geomState(pbDesc.createGeomState()),
  stepState(pbDesc.createGeomState()),
  refState(pbDesc.createGeomState()),
  velocity(pbDesc.solVecInfo()),
  acceleration(pbDesc.solVecInfo()),
  inc_displac(pbDesc.solVecInfo()),
  gravityForce(pbDesc.solVecInfo()),
  elementInternalForce(pbDesc.elemVecInfo()),
  residual(pbDesc.solVecInfo()),
  rhs(pbDesc.solVecInfo()),
  aeroForce(pbDesc.solVecInfo()),
  prev_int_force(pbDesc.solVecInfo()),
  external_force(pbDesc.solVecInfo()),
  stateIncr(pbDesc.solVecInfo()),
  dummyVp(pbDesc.solVecInfo()),
  localDt(pbDesc.getDt()),
  localDelta(pbDesc.getDelta()),
  numStages(pbDesc.getNumStages()),
  maxNumIter(pbDesc.getMaxit()),
  dlambda(1.0 / pbDesc.getNumStages())
{
  probDesc.getConstForce(gravityForce);
  VecType initialDisplacement(pbDesc.solVecInfo());
  /*int aeroAlg =*/ probDesc.getInitState(initialDisplacement, velocity, acceleration, dummyVp);
  setCurrentDisplacement(initialDisplacement);
  double initialTime;
  probDesc.getInitialTime(currStep, initialTime);
  currentTimeIs(initialTime);  
  probDesc.getNewmarkParameters(beta, gamma, alphaf, alpham);
}

NLDynamTimeIntegrator::~NLDynamTimeIntegrator()
{
  delete refState;
  delete stepState;
  delete geomState;
}

// Public methods
// --------------
void
NLDynamTimeIntegrator::currentTimeIs(double newTime)
{
  currTime = newTime;
  midTime = newTime + localDt;
  rhs.zero();
  aeroForce.zero();
  updateForce();
}

// Private methods
// -----------------
void
NLDynamTimeIntegrator::updateForce()
{
  probDesc.getExternalForce(external_force, gravityForce, -1, currTime, geomState, elementInternalForce, aeroForce);
  residual = external_force;
  probDesc.getStiffAndForce(*geomState, residual, elementInternalForce, midTime);
  prev_int_force = external_force - residual;
}

// Time integration
// ----------------
void NLDynamTimeIntegrator::integrate(int numSteps)
{
  double currentRes, resN;
  probDesc.reBuild(*geomState, 0, localDelta); 
  // Begin time-marching
  for (int lastStep = currStep + numSteps; currStep < lastStep; ++currStep)
  {
    currTime += localDt;
    midTime = currTime - localDelta;
    if (verboseFlag)
    {
       fprintf(stderr,"\n**** Begin Time Step %d - Time %e ****\n", currStep + 1, currTime);
    }
    probDesc.getExternalForce(external_force, gravityForce, 1 /* "currStep" */, midTime, geomState, elementInternalForce, aeroForce);
    residual = external_force - prev_int_force;
    stateIncr.zero();
    probDesc.formRHSpredictor(velocity, acceleration, residual, rhs, *geomState, midTime, localDelta); // rhs = [M] (disp + delta * velocity) + delta^2 * residual
    resN = rhs.norm();
    currentRes = resN;
    residual = rhs;
    probDesc.getSolver()->reSolve(rhs); // Solve ([M] + delta^2 [Kt]) du = rhs
    int converged = probDesc.checkConvergence(0, resN, residual, rhs, currTime);
    double lambda = 0.0;
    for (int iStage = 0; iStage < numStages; ++iStage)
    {
      lambda += dlambda;

      for (int iter = 0; iter < maxNumIter; ++iter)
      {
        residual = lambda * external_force;
        geomState->update(stateIncr);
        probDesc.getStiffAndForce(*geomState, residual, elementInternalForce, midTime); // Update elementary matrices & residual force
        prev_int_force = lambda * external_force - residual;
        probDesc.reBuild(*geomState, iter + 1, localDelta); // Assemble [Kt] and factor ([M] + delta^2 * [Kt])
        geomState->get_inc_displacement(inc_displac, *stepState); // Compute incremental displacement
        resN = probDesc.formRHScorrector(inc_displac, velocity, acceleration, residual, rhs, localDelta); // rhs = delta^2 * residual - [M] (inc_displac - delta * velocity) 
        if (verboseFlag) fprintf(stderr,"2 NORMS: fext*fext %e residual*residual %e\n", external_force * external_force, resN * resN);
        currentRes = resN;
        residual = rhs;
        probDesc.getSolver()->reSolve(rhs); // Solve ([M] + delta^2 [Kt]) du = rhs
        converged = probDesc.checkConvergence(iter + 1, resN, residual, rhs, currTime);
        stateIncr = rhs;
        if (converged == 1)
        {
          break;
        }
        else if (converged == -1)
        {
          // Dirty exit...
          fprintf(stderr," ... Exiting, Solution diverged\n");
          exit(-1);
        }
      }
      if (numStages != 1)
      {
        fprintf(stderr,"======> End Stage %d\n", iStage + 1);
      }
    }
    if (converged == 0)
    {
      fprintf(stderr," *** WARNING: Newton solver did not reach convergence after %d iterations (res = %e, target = %e)\n", maxNumIter, currentRes, probDesc.getTolerance());
    }
    geomState->midpoint_step_update(velocity, acceleration, localDelta, *stepState, beta, gamma, alphaf, alpham);
    // if (step+1 == maxStep)  probDesc->processLastOutput(); // Was I right to deactivate ?
    postProcessor().dynamOutput(geomState, velocity, dummyVp, currTime, currStep, external_force, aeroForce, acceleration);
  }
}
