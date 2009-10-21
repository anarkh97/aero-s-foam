#include "JumpProjectorImpl.h"

namespace Pita {

JumpProjectorImpl::JumpProjectorImpl(const String & name,
                                     const JumpProjectorImpl::Manager * manager) :
  JumpProjector(name),
  manager_(manager)
{}

size_t
JumpProjectorImpl::reducedBasisSize() const {
  return reducedBasis()->stateCount();
}

void
JumpProjectorImpl::predictedSeedIs(const Seed * ps) {
  setPredictedSeed(ps);
}

void
JumpProjectorImpl::actualSeedIs(const Seed * as) {
  setActualSeed(as);
}

void
JumpProjectorImpl::seedJumpIs(Seed * sj) {
  setSeedJump(sj);
}

void
JumpProjectorImpl::iterationIs(IterationRank i) {
  DynamState newJumpSeed = actualSeed()->state() - predictedSeed()->state();
  
  seedJump()->statusIs(actualSeed()->status());
  seedJump()->stateIs(newJumpSeed);
  seedJump()->iterationIs(actualSeed()->iteration());

  Vector projectionResult(reducedBasisSize());

  //log() << "seedJump()->state().vectorSize() = " << seedJump()->state().vectorSize() << "\n";

  for (size_t index = 0; index < reducedBasisSize(); ++index) {
    //log() << "reducedBasis()->state(index).vectorSize() = " << reducedBasis()->state(index).vectorSize() << "\n";
    projectionResult[index] = seedJump()->state() * reducedBasis()->state(index); 
  }

  reducedSeedJump()->stateIs(projectionResult);
  reducedSeedJump()->iterationIs(seedJump()->iteration());
  reducedSeedJump()->statusIs(seedJump()->status());
  setIteration(i); 
}

JumpProjectorImpl::Manager::Manager(const DynamStateBasis * drb) :
  defaultReducedBasis_(drb)
{}

JumpProjectorImpl * 
JumpProjectorImpl::Manager::createNewInstance(const String & key) {
  String instanceName = String("JumpProjector ") + key;
  return new JumpProjectorImpl(instanceName, this);
}

} // end namespace Pita
