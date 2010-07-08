#include "ConcurrentBasisManager.h"

namespace Pita {

namespace ConcurrentBasis {

Manager::Manager() :
  linearizedPropagator_(LinPropFactory())
{}

const Propagator *
Manager::concurrentPropagator(const DynamStatePlainBasis * concurrentBasis) const {
   ReactorMap::const_iterator it = reactor_.find(concurrentBasis);
   return (it != reactor_.end()) ? static_cast<const Propagator *>(it->second->notifier()) : NULL;
}

void
Manager::concurrentPropagatorIs(DynamStatePlainBasis * concurrentBasis, const Propagator * cp) {
  if (!concurrentBasis) return;
  if (!cp) {
    reactor_.erase(concurrentBasis);
  } else {
    reactor_[concurrentBasis] = new PropagatorReactor(cp, concurrentBasis, this);
  }
}

PropagatorReactor::PropagatorReactor(const Propagator * notifier, DynamStatePlainBasis * basis, Manager * parent) :
  DynamPropagator::Notifiee(notifier),
  parent_(parent),
  basis_(basis),
  stepReactor_(NULL)
{}

void
PropagatorReactor::onInitialState() {
  const Propagator * actualNotifier = static_cast<const Propagator *>(notifier()); // Safe cast
  NlDynamTimeIntegrator::PtrConst integrator = actualNotifier->integrator();
  
  LinearizedPropagator::Ptr stepPropagator = parent_->linearizedPropagator_.instance(integrator.ptr());
  if (!stepPropagator) {
    stepPropagator = parent_->linearizedPropagator_.instanceNew(integrator.ptr());
  }
   
  stepReactor_ = new IntegratorReactor(integrator.ptr(), basis_.ptr(), stepPropagator.ptr()); 
}

void
PropagatorReactor::onFinalState() {
  stepReactor_ = NULL;
}

IntegratorReactor::IntegratorReactor(const NlDynamTimeIntegrator * notifier,
                                     DynamStatePlainBasis * basis,
                                     DynamPropagator * stepPropagator) :
  DynamTimeIntegrator::NotifieeConst(notifier),
  basis_(basis),
  stepPropagator_(stepPropagator)
{}

void
IntegratorReactor::onCurrentCondition() {
  size_t stateCount = basis_->stateCount();
  for (size_t i = 0; i < stateCount; ++i) {
    stepPropagator_->initialStateIs(basis_->state(i));
    basis_->stateIs(i, stepPropagator_->finalState());
  }
}

} /* end namespace ConcurrentBasis */

} /* end namespace Pita */
