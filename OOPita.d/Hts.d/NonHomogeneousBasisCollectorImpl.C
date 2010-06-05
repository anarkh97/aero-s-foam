#include "NonHomogeneousBasisCollectorImpl.h"

namespace Pita { namespace Hts {

NonHomogeneousBasisCollectorImpl::NonHomogeneousBasisCollectorImpl() :
  BasisCollectorImpl()
{}

DynamState
NonHomogeneousBasisCollectorImpl::nonHomogeneousState(const HalfSliceId & sliceId) const {
  PropagationReactorContainer::const_iterator it = this->propagationReactor().find(sliceId);
  if (it != this->propagationReactor().end() && it->first == sliceId) {
    // TODO Remove dynamic_cast
    const PropagationReactor * r = dynamic_cast<PropagationReactor *>(it->second.ptr());
    return r->nonHomogeneousState();
  } else {
    return DynamState();
  }
}

BasisCollectorImpl::PropagationReactor *
NonHomogeneousBasisCollectorImpl::propagationReactorNew(
    const DynamPropagator * notifier,
    const HalfSliceId & id) {
  return new PropagationReactor(notifier, id, this);
}

NonHomogeneousBasisCollectorImpl::PropagationReactor::PropagationReactor(
    const DynamPropagator * notifier,
    const HalfSliceId & id,
    NonHomogeneousBasisCollectorImpl * parent) :
  BasisCollectorImpl::PropagationReactor(notifier, id, parent),
  nonHomogeneousState_()
{}

void
NonHomogeneousBasisCollectorImpl::PropagationReactor::onFinalState() {
  if (nonHomogeneousState().vectorSize() != 0) {
    parent()->finalStateIs(sliceId(), notifier()->finalState());
  } else {
    nonHomogeneousState_ = notifier()->finalState();
  }
}

} /* end namespace Pita */ } /* end namespace Hts */
