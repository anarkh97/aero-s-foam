#include "BasisCollectorImpl.h"

namespace Pita { namespace Hts {

BasisCollectorImpl::BasisCollectorImpl() :
  BasisCollector(),
  propagationReactor_()
{}

IntegratorPropagator *
BasisCollectorImpl::source(const HalfSliceId & sliceId) const {
  PropagationReactorContainer::const_iterator it = propagationReactor_.find(sliceId);
  return (it != propagationReactor_.end()) ?
    static_cast<IntegratorPropagator *>(const_cast<DynamPropagator *>(it->second->notifier())) :
    NULL;
}

size_t
BasisCollectorImpl::sourceCount() const {
  return propagationReactor_.size();
}

IntegratorPropagator *
BasisCollectorImpl::sourceNew(const HalfSliceId & sliceId, DynamTimeIntegrator * baseIntegrator) {
  // Find insertion point
  PropagationReactorContainer::iterator it = propagationReactor_.lower_bound(sliceId);
  if (it != propagationReactor_.end() && it->first == sliceId)
    throw NameInUseException();
 
  // Build propagator and add new reactor
  IntegratorPropagator::Ptr newPropagator = IntegratorPropagator::New(baseIntegrator);
  PropagationReactor::Ptr newReactor = propagationReactorNew(newPropagator.ptr(), sliceId);

  propagationReactor_.insert(it, std::make_pair(sliceId, newReactor));
    
  return newPropagator.ptr();
}

void
BasisCollectorImpl::sourceDel(const HalfSliceId & sliceId) {
  propagationReactor_.erase(sliceId);
}

BasisCollectorImpl::CollectedState
BasisCollectorImpl::firstForwardFinalState() const {
  return !forwardFinalState_.empty() ? forwardFinalState_.top() : make_pair(HalfSliceRank(-1), DynamState());
}

void
BasisCollectorImpl::firstForwardFinalStateDel() {
  forwardFinalState_.pop();
}

BasisCollectorImpl::CollectedState
BasisCollectorImpl::firstBackwardFinalState() const {
  return !backwardFinalState_.empty() ? backwardFinalState_.top() : make_pair(HalfSliceRank(-1), DynamState());
}

void
BasisCollectorImpl::firstBackwardFinalStateDel() {
  backwardFinalState_.pop();
}

void
BasisCollectorImpl::finalStateIs(const HalfSliceId & sliceId, const DynamState & state) {
  HalfTimeSlice::Direction dir = sliceId.direction();
  StateContainer & stateContainer = (dir == HalfTimeSlice::FORWARD) ? forwardFinalState_ : backwardFinalState_;
  stateContainer.push(make_pair(sliceId.rank(), state));
}

BasisCollectorImpl::PropagationReactor *
BasisCollectorImpl::propagationReactorNew(IntegratorPropagator * notifier,
                                                   const HalfSliceId & id) {
  return new PropagationReactor(notifier, id, this);
}

BasisCollectorImpl::PropagationReactor::PropagationReactor(DynamPropagator * notifier,
                                                                    const HalfSliceId & id,
                                                                    BasisCollectorImpl * parent) :
  DynamPropagator::Notifiee(notifier),
  sliceId_(id),
  parent_(parent)
{}

void
BasisCollectorImpl::PropagationReactor::onFinalState() {
  parent_->finalStateIs(sliceId_, notifier()->finalState());
}

} /* end notifier Hts */ } /* end namespace Pita */
