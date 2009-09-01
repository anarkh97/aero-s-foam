#include "JumpProjectorImpl.h"

namespace Pita {

JumpProjectorImpl::JumpProjectorImpl(const JumpProjectorImpl::Manager * manager,
                                     JumpBuilder * jumpBuilder) :
  manager_(manager),
  jumpBuilder_(jumpBuilder)
{}

size_t
JumpProjectorImpl::reducedBasisSize() const {
  return reducedBasis()->stateCount();
}

void
JumpProjectorImpl::predictedSeedIs(const Seed * ps) {
  //log() << "JP: new predicted seed " << ps->name() << "\n";
  jumpBuilder_->seedIs(JumpBuilder::RIGHT, ps);
  setPredictedSeed(ps);
}

void
JumpProjectorImpl::actualSeedIs(const Seed * as) {
  //log() << "JP: new actual seed " << as->name() << "\n";
  jumpBuilder_->seedIs(JumpBuilder::LEFT, as);
  setActualSeed(as);
}

void
JumpProjectorImpl::seedJumpIs(Seed * sj) {
  //log() << "JP: new jump seed " << sj->name() << "\n";
  jumpBuilder_->jumpSeedIs(sj);
  setSeedJump(sj);
}

void
JumpProjectorImpl::iterationIs(IterationRank i) {
  Vector projectionResult(reducedBasisSize());

  //log() << "seedJump()->state().vectorSize() = " << seedJump()->state().vectorSize() << "\n";

  // TODO optimize
  for (size_t index = 0; index < reducedBasisSize(); ++index) {
    //log() << "reducedBasis()->state(index).vectorSize() = " << reducedBasis()->state(index).vectorSize() << "\n";
    projectionResult[index] = seedJump()->state() * reducedBasis()->state(index); 
  }

  reducedSeedJump()->stateIs(projectionResult);
  reducedSeedJump()->iterationIs(seedJump()->iteration());
  setIteration(i); 
}

JumpProjectorImpl::Manager::Manager(const DynamStateBasis * drb) :
  defaultReducedBasis_(drb)
{}

JumpProjectorImpl * 
JumpProjectorImpl::Manager::createNewInstance(const String & key) {
  JumpBuilder::Ptr jumpBuilder = JumpBuilder::New();
  return new JumpProjectorImpl(this, jumpBuilder.ptr());
}

} // end namespace Pita
