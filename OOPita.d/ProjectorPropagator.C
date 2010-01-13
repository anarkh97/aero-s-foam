#include "ProjectorPropagator.h"

#include "DynamStateOps.h"

#include <cmath>

namespace Pita {

const double ProjectorPropagator::tolerance_ = 1.0e-6;
  
ProjectorPropagator::ProjectorPropagator(DynamPropagator * basePropagator, NlDynamOps * dynamOps) :
  DynamPropagator(basePropagator->vectorSize()),
  basePropagator_(basePropagator),
  dynamOps_(dynamOps),
  refDisp_(basePropagator->vectorSize(), 0.0),
  refTime_(Seconds(0.0)),
  projectionBasis_(DynamStatePlainBasis::New(basePropagator->vectorSize())),
  orthoBasis_(DynamStatePlainBasis::New(basePropagator->vectorSize())),
  propagatedBasis_(DynamStatePlainBasis::New(basePropagator->vectorSize()))
{}

void
ProjectorPropagator::initialStateIs(const DynamState & initialState) {
  DynamState finalState(initialState.vectorSize(), 0.0);
  log() << "Projection on subspace of size " << orthoBasis_->stateCount() << "\n";
  size_t states = orthoBasis_->stateCount();
  for (size_t i = 0; i < states; ++i) {
    double coef = orthoBasis_->state(i) * initialState;
    log() << coef << " ";
    finalState.linAdd(coef, propagatedBasis_->state(i));
  }
  log() << "\n";

  // Commit changes
  setInitialState(initialState);
  initialStateNotify();
  setFinalState(finalState);
  finalStateNotify();
}

void
ProjectorPropagator::projectionBasisIs(DynamStateBasis::PtrConst basis) {
  orthoBasis_->stateBasisDel();
  propagatedBasis_->stateBasisDel();
  projectionBasis_->stateBasisDel();
  projectionBasisInc(basis);
}

void
ProjectorPropagator::referenceDisplacementIs(const GenVector<double> & disp, Seconds time) {
  refDisp_ = disp;
  refTime_ = time;
  if (projectionBasis_->stateCount() > 0) {
    DynamStatePlainBasis::Ptr tempBasis = projectionBasis_;
    projectionBasis_ = DynamStatePlainBasis::New(tempBasis->vectorSize());
    orthoBasis_->stateBasisDel();
    propagatedBasis_->stateBasisDel();
    projectionBasisInc(tempBasis);
  }
}

void
ProjectorPropagator::projectionBasisInc(DynamStateBasis::PtrConst lastBasis) {
  log() << "ProjectorPropagator::projectionBasisInc - Current = " << projectionBasis_->stateCount() << " - Adding = " << lastBasis->stateCount();
  if (lastBasis->stateCount() != 0) {

    dynamOps_->displacementIs(refDisp_, refTime_);

    size_t statesToAdd = lastBasis->stateCount();
    for (size_t i = 0; i < statesToAdd; ++i) {
      DynamState currentState = lastBasis->state(i);

      size_t statesInBasis = orthoBasis_->stateCount();
      for (size_t j = 0; j < statesInBasis; ++j) {
        double coef = -(currentState * orthoBasis_->state(j));
        currentState.linAdd(coef, projectionBasis_->state(j));
      }

      DynamState orthoState = mult(dynamOps_.ptr(), currentState);

      double norm_i = std::sqrt(currentState * orthoState);
      if (norm_i < tolerance_)
        continue;
      norm_i = 1.0 / norm_i;
      currentState.displacement() *= norm_i;
      currentState.velocity() *= norm_i;
      projectionBasis_->lastStateIs(currentState);
      orthoState.displacement() *= norm_i;
      orthoState.velocity() *= norm_i;
      orthoBasis_->lastStateIs(orthoState);
    } 
  }
  log() << " - Final = " << projectionBasis_->stateCount() << "\n";
  propagatedBasis_->stateBasisDel();
  propagatedBasis_->lastStateBasisIs(projectionBasis_.ptr());
  /*for (int i = 0; i < orthoBasis_->stateCount(); ++i) {
    for (int j = 0; j < propagatedBasis_->stateCount(); ++j) {
       log() << orthoBasis_->state(i) * propagatedBasis_->state(j) << " ";
    }
    log() << "\n";
  }*/
  timeStepCount_ = TimeStepCount(0);
}

void
ProjectorPropagator::timeStepCountInc() {
  size_t stateCount = propagatedBasis_->stateCount();
  for (size_t i = 0; i < stateCount; ++i) {
    basePropagator_->initialStateIs(propagatedBasis_->state(i));
    propagatedBasis_->stateIs(i, basePropagator_->finalState());
  }
  ++timeStepCount_;
}
  
void
ProjectorPropagator::toString(OStream & out) const {
  int orthoBasisSize = orthoBasis_->stateCount();
  int projectionBasisSize = projectionBasis_->stateCount();
  for (int i = 0; i < orthoBasisSize; ++i) {
    for (int j = 0; j < projectionBasisSize; ++j) {
       out << orthoBasis_->state(i) * projectionBasis_->state(j) << " ";
    }
    out << "\n";
  }
}
  
} // end namespace Pita
