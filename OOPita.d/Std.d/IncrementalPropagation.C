#include "IncrementalPropagation.h"

namespace Pita { namespace Std {

IncrementalPropagation::IncrementalPropagation(const Fwk::String & name,
                                               AffineDynamPropagator * propagator) :
  NamedTask(name),
  propagator_(propagator),
  previousSeedState_()
{}

void
IncrementalPropagation::seedIs(const Seed * s) {
  if (s->status() == Seed::ACTIVE) {
    previousSeedState_ = s->state();
  }
  seed_ = s;
}

void
IncrementalPropagation::propagatedSeedIs(Seed * ps) {
  propagatedSeed_ = ps;
}

void
IncrementalPropagation::iterationIs(IterationRank i) {
  // Sanity checks
  assert(seed()->iteration() == i);
  assert(seed()->status() != Seed::INACTIVE);

  // Perform propagation
  if (previousSeedState_.vectorSize() != 0) {
    //log() << "Reuse state\n";
    propagator()->constantTermIs(AffineDynamPropagator::HOMOGENEOUS);
    //log() << "External force flagged off\n"; 
    propagator()->initialStateIs(seed()->state() - previousSeedState_);
  } else {
    //log() << "New state\n";
    propagator()->initialStateIs(seed()->state());
  }
 
  propagatedSeed()->statusIs(seed()->status());
  if (propagatedSeed()->state().vectorSize() != 0) {
    propagatedSeed()->stateIs(propagatedSeed()->state() + propagator()->finalState());
  } else {
    propagatedSeed()->stateIs(propagator()->finalState());
  }

  previousSeedState_ = seed()->state(); 

  propagatedSeed()->iterationIs(seed()->iteration().next());
}

} /* end namespace Std */ } /* end namespace Pita */
