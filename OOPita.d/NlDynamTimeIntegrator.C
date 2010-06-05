#include "NlDynamTimeIntegrator.h"
#include "PitaNonLinDynam.h"
#include <Corotational.d/GeomState.h>
#include <cstdio>
#include <cstdlib>
#include <cmath>

using std::fprintf;
using std::exit;

namespace Pita {

NlDynamTimeIntegrator::NlDynamTimeIntegrator(PitaNonLinDynamic * pbDesc) :
  DynamTimeIntegrator(pbDesc->solVecInfo()),
  probDesc_(pbDesc),
  geomState_(pbDesc->createGeomState()),
  stepState_(pbDesc->createGeomState()),
  refState_(pbDesc->createGeomState()),
  displacement_(pbDesc->solVecInfo()),
  velocity_(pbDesc->solVecInfo()),
  incDisplac_(pbDesc->solVecInfo()),
  gravityForce_(pbDesc->solVecInfo()),
  elementInternalForce_(pbDesc->elemVecInfo()),
  residual_(pbDesc->solVecInfo()),
  rhs_(pbDesc->solVecInfo()),
  aeroForce_(pbDesc->solVecInfo()),
  prevIntForce_(pbDesc->solVecInfo()),
  externalForce_(pbDesc->solVecInfo()),
  stateIncr_(pbDesc->solVecInfo()),
  dummyVp_(pbDesc->solVecInfo()),
  localDt_(pbDesc->getDt()),
  localDelta_(pbDesc->getDelta()),
  numStages_(pbDesc->getNumStages()),
  maxNumIter_(pbDesc->getMaxit()),
  dlambda_(1.0 / pbDesc->getNumStages())
{
  pbDesc->getConstForce(gravityForce_);
  setTimeStepSize(Seconds(localDt_));
  GenVector<double> dummyAcceleration(pbDesc->solVecInfo());
  probDesc_->getInitState(displacement_, velocity_, dummyAcceleration, dummyVp_);
  double initTime;
  probDesc_->getInitialTime(currStep_, initTime);
  initialConditionIs(DynamState(displacement_, velocity_), Seconds(initTime));
  if (currStep_ != 0)
    setTimeStepCount(TimeStepCount(currStep_));
}

NlDynamTimeIntegrator::~NlDynamTimeIntegrator() {
  delete refState_;
  delete stepState_;
  delete geomState_;
}

void
NlDynamTimeIntegrator::timeStepSizeIs(Seconds dt) {
  updateDelta(dt.value());
  setTimeStepSize(dt);
}

void
NlDynamTimeIntegrator::initialConditionIs(const DynamState & initialState, Seconds initialTime) {
  currTime_ = initialTime.value();
  midTime_ = currTime_ + localDelta_;
  rhs_.zero();
  aeroForce_.zero();
  *geomState_ = *refState_;
  geomState_->update(initialState.displacement());
  *stepState_ = *geomState_;
  velocity_ = initialState.velocity();
  updateForce();
  setInitialTime(initialTime);
  setInitialState(initialState);
  setTimeStepCount(TimeStepCount(0));
  setCurrentTime(initialTime);
  setCurrentState(initialState);
  performNotification(&NotifieeConst::onInitialCondition);
}

void
NlDynamTimeIntegrator::currentTimeInc(Seconds increment) {
  unsigned int stepCount = static_cast<unsigned int>(std::floor(increment.value() / localDt_));
  double remainder = increment.value() - stepCount * localDt_;
  integrate(stepCount);
  updateDelta(remainder);
  integrate(1);
  updateDelta(timeStepSize().value());
}

void
NlDynamTimeIntegrator::timeStepCountInc(TimeStepCount steps) {
  integrate(steps.value());
}


// Private methods

void
NlDynamTimeIntegrator::updateDelta(double dt) {
  localDt_ = dt;
  localDelta_ = 0.5 * dt;
  midTime_ = currTime_ + localDelta_;
}

void
NlDynamTimeIntegrator::updateForce() {
  probDesc_->getExternalForce(externalForce_, gravityForce_, -1, currTime_, geomState_, elementInternalForce_, aeroForce_);
  residual_ = externalForce_;
  probDesc_->getStiffAndForce(*geomState_, residual_, elementInternalForce_, midTime_);
  prevIntForce_ = externalForce_ - residual_;
}

void
NlDynamTimeIntegrator::integrate(unsigned int steps) {
  double currentRes, resN;
  probDesc_->reBuild(*geomState_, 0, localDelta_); 
  // Begin time-marching
  for (int lastStep = currStep_ + steps; currStep_ < lastStep; ++currStep_)
  {
    currTime_ += localDt_;
    midTime_ = currTime_ - localDelta_;
    if (verboseFlag)
    {
       fprintf(stderr,"\n**** Begin Time Step %d - Time %e ****\n", currStep_ + 1, currTime_);
    }
    probDesc_->getExternalForce(externalForce_, gravityForce_, 1 /* "currStep" */, midTime_, geomState_, elementInternalForce_, aeroForce_);
    residual_ = externalForce_ - prevIntForce_;
    stateIncr_.zero();
    probDesc_->formRHSpredictor(velocity_, residual_, rhs_, *geomState_, midTime_, localDelta_); // rhs = [M] (disp + delta * velocity) + delta^2 * residual
    resN = rhs_.norm();
    currentRes = resN;
    residual_ = rhs_;
    probDesc_->getSolver()->reSolve(rhs_); // Solve ([M] + delta^2 [Kt]) du = rhs
    int converged = probDesc_->checkConvergence(0, resN, residual_, rhs_, currTime_);
    double lambda = 0.0;
    for (int iStage = 0; iStage < numStages_; ++iStage)
    {
      lambda += dlambda_;

      for (int iter = 0; iter < maxNumIter_; ++iter)
      {
        residual_ = lambda * externalForce_;
        geomState_->update(stateIncr_);
        probDesc_->getStiffAndForce(*geomState_, residual_, elementInternalForce_, midTime_); // Update elementary matrices & residual force
        prevIntForce_ = lambda * externalForce_ - residual_;
        probDesc_->reBuild(*geomState_, iter + 1, localDelta_); // Assemble [Kt] and factor ([M] + delta^2 * [Kt])
        geomState_->get_inc_displacement(incDisplac_, *stepState_); // Compute incremental displacement
        resN = probDesc_->formRHScorrector(incDisplac_, velocity_, residual_, rhs_, localDelta_); // rhs = delta^2 * residual - [M] (inc_displac - delta * velocity) 
        if (verboseFlag)
          fprintf(stderr,"2 NORMS: fext*fext %e residual*residual %e\n", externalForce_ * externalForce_, resN * resN);
        currentRes = resN;
        residual_ = rhs_;
        probDesc_->getSolver()->reSolve(rhs_); // Solve ([M] + delta^2 [Kt]) du = rhs
        converged = probDesc_->checkConvergence(iter + 1, resN, residual_, rhs_, currTime_);
        stateIncr_ = rhs_;
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
      if (numStages_ != 1)
      {
        fprintf(stderr,"======> End Stage %d\n", iStage + 1);
      }
    }
    if (converged == 0)
    {
      fprintf(stderr," *** WARNING: Newton solver did not reach convergence after %d iterations (res = %e, target = %e)\n", maxNumIter_, currentRes, probDesc_->getTolerance());
    }
    geomState_->midpoint_step_update(velocity_, localDelta_, *stepState_);
    // if (step+1 == maxStep)  probDesc_->processLastOutput(); // Was I right to deactivate ?
    
    // Update attributes 
    setCurrentTime(Seconds(currTime_));
    geomState_->get_inc_displacement(displacement_, *refState_, false);
    setCurrentState(DynamState(displacement_, velocity_));
    setTimeStepCount(TimeStepCount(currStep_ + 1));
    // Perform notification (To allow side-effects, such as post-processing)
    performNotification(&NotifieeConst::onCurrentCondition);
  }
}

} // end namespace Pita
