#include "HalfSliceBasisCollectorImpl.h"

namespace Pita {

HalfSliceBasisCollectorImpl::HalfSliceBasisCollectorImpl() :
  HalfSliceBasisCollector(),
  propagationReactor_()
{}

IntegratorPropagator *
HalfSliceBasisCollectorImpl::source(const HalfSliceId & sliceId) const {
  PropagationReactorContainer::const_iterator it = propagationReactor_.find(sliceId);
  return (it != propagationReactor_.end()) ?
    static_cast<IntegratorPropagator *>(const_cast<DynamPropagator *>(it->second->notifier())) :
    NULL;
}

size_t
HalfSliceBasisCollectorImpl::sourceCount() const {
  return propagationReactor_.size();
}

IntegratorPropagator *
HalfSliceBasisCollectorImpl::sourceNew(const HalfSliceId & sliceId, DynamTimeIntegrator * baseIntegrator) {
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
HalfSliceBasisCollectorImpl::sourceDel(const HalfSliceId & sliceId) {
  propagationReactor_.erase(sliceId);
}

HalfSliceBasisCollectorImpl::CollectedState
HalfSliceBasisCollectorImpl::firstForwardFinalState() const {
  return !forwardFinalState_.empty() ? forwardFinalState_.top() : make_pair(HalfSliceRank(-1), DynamState());
}

void
HalfSliceBasisCollectorImpl::firstForwardFinalStateDel() {
  forwardFinalState_.pop();
}

HalfSliceBasisCollectorImpl::CollectedState
HalfSliceBasisCollectorImpl::firstBackwardFinalState() const {
  return !backwardFinalState_.empty() ? backwardFinalState_.top() : make_pair(HalfSliceRank(-1), DynamState());
}

void
HalfSliceBasisCollectorImpl::firstBackwardFinalStateDel() {
  backwardFinalState_.pop();
}

void
HalfSliceBasisCollectorImpl::finalStateIs(const HalfSliceId & sliceId, const DynamState & state) {
  HalfTimeSlice::Direction dir = sliceId.direction();
  StateContainer & stateContainer = (dir == HalfTimeSlice::FORWARD) ? forwardFinalState_ : backwardFinalState_;
  stateContainer.push(make_pair(sliceId.rank(), state));
}

HalfSliceBasisCollectorImpl::PropagationReactor *
HalfSliceBasisCollectorImpl::propagationReactorNew(IntegratorPropagator * notifier,
                                                   const HalfSliceId & id) {
  return new PropagationReactor(notifier, id, this);
}

HalfSliceBasisCollectorImpl::PropagationReactor::PropagationReactor(DynamPropagator * notifier,
                                                                    const HalfSliceId & id,
                                                                    HalfSliceBasisCollectorImpl * parent) :
  DynamPropagator::Notifiee(notifier),
  sliceId_(id),
  parent_(parent)
{}

void
HalfSliceBasisCollectorImpl::PropagationReactor::onFinalState() {
  parent_->finalStateIs(sliceId_, notifier()->finalState());
}

} // end namespace Pita
