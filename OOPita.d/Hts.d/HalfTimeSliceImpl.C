#include "HalfTimeSliceImpl.h"

#include "../IntegratorPropagator.h"

// TODO Remove ugly HACK
#include "../IntegratorPropagator.h"
#include "../HomogeneousGenAlphaIntegrator.h"

#include <memory>

#include <cassert>

namespace Pita { namespace Hts {

/* HalfTimeSliceImpl */

HalfTimeSliceImpl::HalfTimeSliceImpl(HalfSliceRank r,
                                     HalfTimeSlice::Direction d,
                                     DynamPropagator * dp) :
  HalfTimeSlice(r, d),
  propagator_(dp),
  previousSeedState_()
{}

void
HalfTimeSliceImpl::seedIs(const Seed * s) {
  previousSeedState_ = s->state();
  setSeed(s);
}

void
HalfTimeSliceImpl::iterationIs(IterationRank i) {
  assert(seed()->iteration() == i);
  propagateSeed();
  propagatedSeed()->iterationIs(seed()->iteration().next());
}

void
HalfTimeSliceImpl::propagatedSeedIs(Seed * ps) {
  setPropagatedSeed(ps);
}

void
HalfTimeSliceImpl::propagateSeed() {
  if (seed()) {
    if (previousSeedState_.vectorSize() != 0) {
      //log() << "Reuse state\n";

      // TODO Remove ugly HACK
      if (IntegratorPropagator * prop = dynamic_cast<IntegratorPropagator *>(propagator())) {
        if (AffineGenAlphaIntegrator * integr = dynamic_cast<AffineGenAlphaIntegrator *>(prop->integrator())) {
          integr->externalForceFlagIs(false);
          //log() << "External force flagged off\n"; 
        }
      }
      propagator()->initialStateIs(seed()->state() - previousSeedState_);
    } else {
      //log() << "New state\n";
      propagator()->initialStateIs(seed()->state());
    }
    if (propagatedSeed()) {
      propagatedSeed()->statusIs(seed()->status());
      if (propagatedSeed()->state().vectorSize() != 0) {
        propagatedSeed()->stateIs(propagatedSeed()->state() + propagator()->finalState());
      } else {
        propagatedSeed()->stateIs(propagator()->finalState());
      }
    }
    
    previousSeedState_ = seed()->state(); 
  }
}


/* HalfTimeSliceImpl::Manager */

HalfTimeSliceImpl::Manager::Manager(HalfTimeSliceImpl::Manager::PropagatorManager * propagatorManager) :
  propagatorManager_(propagatorManager)
{}

HalfTimeSliceImpl *
HalfTimeSliceImpl::Manager::instance(const HalfSliceId & id) const {
  return Impl::instance(id);
}

size_t
HalfTimeSliceImpl::Manager::instanceCount() const {
  return Impl::instanceCount();
}

HalfTimeSliceImpl *
HalfTimeSliceImpl::Manager::instanceNew(const HalfSliceId & id) {
  return Impl::instanceNew(id);
}

void
HalfTimeSliceImpl::Manager::instanceDel(const HalfSliceId & id) {
  Impl::instanceDel(id);
}

HalfTimeSliceImpl *
HalfTimeSliceImpl::Manager::createNewInstance(const HalfSliceId & id) {
  // Instantiate new HalfTimeSlice
  DynamPropagator::Ptr newPropagator = propagatorManager_->instanceNew(id);
  std::auto_ptr<HalfTimeSliceImpl> halfSlice(new HalfTimeSliceImpl(id.rank(), id.direction(), newPropagator.ptr()));
  
  return halfSlice.release();
}

} /* end namespace Hts */ } /* end namespace Pita */
