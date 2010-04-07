#include "LocalCorrectionTimeSlice.h"

#include <cassert>

namespace Pita { namespace Hts {

/* LocalCorrectionTimeSlice implementation */
  
LocalCorrectionTimeSlice::LocalCorrectionTimeSlice(HalfSliceRank headRank, DynamPropagator * propagator) :
  CorrectionTimeSlice(headRank),
  //jumpBuilder_(JumpBuilder::New()),
  propagator_(propagator)
{}

void
LocalCorrectionTimeSlice::predictedSeedIs(const Seed * ps) {
  //jumpBuilder_->seedIs(JumpBuilder::RIGHT, ps);
  setPredictedSeed(ps);
}

void
LocalCorrectionTimeSlice::actualSeedIs(const Seed * as) {
  //jumpBuilder_->seedIs(JumpBuilder::LEFT, as);
  setActualSeed(as);
}

void
LocalCorrectionTimeSlice::jumpIs(Seed * j) {
  //jumpBuilder_->jumpSeedIs(j);
  setJump(j);
}

void
LocalCorrectionTimeSlice::iterationIs(IterationRank ir) {
  assert(correction()->iteration() == ir);
  assert(jump()->iteration() == ir);
  
  DynamState seedUpdate = correction()->state() + jump()->state();
  propagator_->initialStateIs(seedUpdate);
  
  Seed::Status updateStatus = (jump()->status() == Seed::CONVERGED) ? correction()->status() : Seed::ACTIVE;
  nextCorrection()->statusIs(updateStatus);
  nextCorrection()->stateIs(propagator_->finalState());
  nextCorrection()->iterationIs(ir);
 
  setIteration(ir);
}

/* Manager implementation */

LocalCorrectionTimeSlice::Manager::Manager(DynamPropagator * sharedPropagator) :
  sharedPropagator_(sharedPropagator)
{}

LocalCorrectionTimeSlice *
LocalCorrectionTimeSlice::Manager::createNewInstance(const HalfSliceRank & key) {
  return new LocalCorrectionTimeSlice(key, sharedPropagator()); 
}

} /* end namespace Hts */ } /* end namespace Pita */
