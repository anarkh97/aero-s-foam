#include "BasisCollectorImpl.h"

namespace Pita { namespace Hts {

BasisCollectorImpl::BasisCollectorImpl() :
  BasisCollector(),
  propagationReactor_()
{}

const DynamPropagator *
BasisCollectorImpl::source(const HalfSliceId & sliceId) const {
  PropagationReactorContainer::const_iterator it = propagationReactor_.find(sliceId);
  return (it != propagationReactor_.end()) ? it->second->notifier() : NULL;
}

size_t
BasisCollectorImpl::sourceCount() const {
  return propagationReactor_.size();
}

void
BasisCollectorImpl::sourceIs(const HalfSliceId & sliceId, const DynamPropagator * source) {
  if (source == NULL) {
    propagationReactor_.erase(sliceId); // Remove current
    return;
  }

  PropagationReactor::Ptr newReactor = propagationReactorNew(source, sliceId);
  
  // Find insertion point
  PropagationReactorContainer::iterator it = propagationReactor_.lower_bound(sliceId);
  if (it != propagationReactor_.end() && it->first == sliceId) {
    it->second = newReactor; // Replace previous
  } else {
    propagationReactor_.insert(it, std::make_pair(sliceId, newReactor)); // Insert new
  }
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
  Direction dir = sliceId.direction();
  StateContainer & stateContainer = (dir == FORWARD) ? forwardFinalState_ : backwardFinalState_;
  stateContainer.push(make_pair(sliceId.rank(), state));
}

BasisCollectorImpl::PropagationReactor *
BasisCollectorImpl::propagationReactorNew(const DynamPropagator * notifier,
                                          const HalfSliceId & id) {
  return new PropagationReactor(notifier, id, this);
}

BasisCollectorImpl::PropagationReactor::PropagationReactor(const DynamPropagator * notifier,
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
