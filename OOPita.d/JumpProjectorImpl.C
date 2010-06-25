#include "JumpProjectorImpl.h"

#include <cassert>

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
JumpProjectorImpl::iterationIs(IterationRank ir) {
  assert(actualSeed()->iteration() == ir);

  updateJump();
  //log() << "seedJump()->state().vectorSize() = " << seedJump()->state().vectorSize() << "\n";

  Vector projectionResult(reducedBasisSize());
  for (size_t index = 0; index < reducedBasisSize(); ++index) {
    //log() << "reducedBasis()->state(index).vectorSize() = " << reducedBasis()->state(index).vectorSize() << "\n";
    projectionResult[index] = seedJump()->state() * reducedBasis()->state(index); 
  }

  reducedSeedJump()->stateIs(projectionResult);
  reducedSeedJump()->iterationIs(seedJump()->iteration());
  reducedSeedJump()->statusIs(seedJump()->status());

  setIteration(ir); 
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
