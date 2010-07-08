#include "BasicHalfTimeSliceImpl.h"

#include <cassert>

namespace Pita { namespace Hts {

/* BasicHalfTimeSliceImpl */

BasicHalfTimeSliceImpl::BasicHalfTimeSliceImpl(HalfSliceRank r,
                                               HalfTimeSlice::Direction d,
                                               DynamPropagator * dp) :
  HalfTimeSlice(r, d),
  propagator_(dp)
{}

void
BasicHalfTimeSliceImpl::iterationIs(IterationRank i) {
  assert(seed()->iteration() == i);
  assert(seed()->status() != Seed::INACTIVE);
  propagateSeed();
  propagatedSeed()->iterationIs(seed()->iteration().next());
}

void
BasicHalfTimeSliceImpl::propagateSeed() {
  propagator()->initialStateIs(seed()->state());
 
  propagatedSeed()->stateIs(propagator()->finalState());
  propagatedSeed()->statusIs(seed()->status());
}


/* BasicHalfTimeSliceImpl::Manager */

BasicHalfTimeSliceImpl::Manager::Manager(BasicHalfTimeSliceImpl::Manager::DynamPropagatorManager * propagatorManager) :
  propagatorManager_(propagatorManager)
{}

BasicHalfTimeSliceImpl *
BasicHalfTimeSliceImpl::Manager::instance(const HalfSliceId & id) const {
  return Impl::instance(id);
}

size_t
BasicHalfTimeSliceImpl::Manager::instanceCount() const {
  return Impl::instanceCount();
}

BasicHalfTimeSliceImpl *
BasicHalfTimeSliceImpl::Manager::instanceNew(const HalfSliceId & id) {
  return Impl::instanceNew(id);
}

void
BasicHalfTimeSliceImpl::Manager::instanceDel(const HalfSliceId & id) {
  Impl::instanceDel(id);
}

BasicHalfTimeSliceImpl *
BasicHalfTimeSliceImpl::Manager::createNewInstance(const HalfSliceId & id) {
  DynamPropagator::Ptr newPropagator = propagatorManager_->instanceNew(id);
  return new BasicHalfTimeSliceImpl(id.rank(), id.direction(), newPropagator.ptr());
}

} /* end namespace Hts */ } /* end namespace Pita */
