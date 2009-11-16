#include "AffineHalfTimeSliceImpl.h"

#include <memory>
#include <cassert>

namespace Pita { namespace Hts {

/* AffineHalfTimeSliceImpl */

AffineHalfTimeSliceImpl::AffineHalfTimeSliceImpl(HalfSliceRank r,
                                     HalfTimeSlice::Direction d,
                                     AffineDynamPropagator * dp) :
  HalfTimeSlice(r, d),
  propagator_(dp),
  previousSeedState_()
{}

void
AffineHalfTimeSliceImpl::seedIs(const Seed * s) {
  previousSeedState_ = s->state();
  setSeed(s);
}

void
AffineHalfTimeSliceImpl::iterationIs(IterationRank i) {
  assert(seed()->iteration() == i);
  assert(seed()->status() != Seed::INACTIVE);
  propagateSeed();
  propagatedSeed()->iterationIs(seed()->iteration().next());
}

void
AffineHalfTimeSliceImpl::propagatedSeedIs(Seed * ps) {
  setPropagatedSeed(ps);
}

void
AffineHalfTimeSliceImpl::propagateSeed() {
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
}


/* AffineHalfTimeSliceImpl::Manager */

AffineHalfTimeSliceImpl::Manager::Manager(AffineHalfTimeSliceImpl::Manager::PropagatorManager * propagatorManager) :
  propagatorManager_(propagatorManager)
{}

AffineHalfTimeSliceImpl *
AffineHalfTimeSliceImpl::Manager::instance(const HalfSliceId & id) const {
  return Impl::instance(id);
}

size_t
AffineHalfTimeSliceImpl::Manager::instanceCount() const {
  return Impl::instanceCount();
}

AffineHalfTimeSliceImpl *
AffineHalfTimeSliceImpl::Manager::instanceNew(const HalfSliceId & id) {
  return Impl::instanceNew(id);
}

void
AffineHalfTimeSliceImpl::Manager::instanceDel(const HalfSliceId & id) {
  Impl::instanceDel(id);
}

AffineHalfTimeSliceImpl *
AffineHalfTimeSliceImpl::Manager::createNewInstance(const HalfSliceId & id) {
  // Instantiate new HalfTimeSlice
  AffineDynamPropagator::Ptr newPropagator = propagatorManager_->instanceNew(id);
  std::auto_ptr<AffineHalfTimeSliceImpl> halfSlice(new AffineHalfTimeSliceImpl(id.rank(), id.direction(), newPropagator.ptr()));
  
  return halfSlice.release();
}

} /* end namespace Hts */ } /* end namespace Pita */
