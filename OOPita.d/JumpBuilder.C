#include "JumpBuilder.h"

#include <cassert>

namespace Pita {

void
JumpBuilder::iterationIs(IterationRank ir) {
  assert(actualSeed()->iteration().next() == ir);
  
  updateJump();
  setIteration(ir);
} 

void
JumpBuilder::updateJump() {
  assert(actualSeed()->iteration() == predictedSeed()->iteration() || actualSeed()->iteration() == predictedSeed()->iteration().next());
  assert(actualSeed()->status() != Seed::INACTIVE);
  assert(actualSeed()->status() != Seed::SPECIAL);

  DynamState newJumpSeed = actualSeed()->state() - predictedSeed()->state();
  
  seedJump()->statusIs(actualSeed()->status());
  seedJump()->stateIs(newJumpSeed);
  seedJump()->iterationIs(actualSeed()->iteration());
}

JumpBuilder *
JumpBuilder::ManagerImpl::createNewInstance(const String & key) {
  String taskName = "JumpBuilder " + toString(key);
  return new JumpBuilder(taskName);
}

} /* end namespace */
