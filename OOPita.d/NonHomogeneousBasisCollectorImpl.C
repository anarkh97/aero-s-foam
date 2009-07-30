#include "NonHomogeneousBasisCollectorImpl.h"

namespace Pita { namespace Hts {

NonHomogeneousBasisCollectorImpl::NonHomogeneousBasisCollectorImpl() :
  HalfSliceBasisCollectorImpl()
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

HalfSliceBasisCollectorImpl::PropagationReactor *
NonHomogeneousBasisCollectorImpl::propagationReactorNew(
    IntegratorPropagator * notifier,
    const HalfSliceId & id) {
  return new PropagationReactor(notifier, id, this);
}

NonHomogeneousBasisCollectorImpl::PropagationReactor::PropagationReactor(
    IntegratorPropagator * notifier,
    const HalfSliceId & id,
    NonHomogeneousBasisCollectorImpl * parent) :
  HalfSliceBasisCollectorImpl::PropagationReactor(notifier, id, parent),
  nonHomogeneousState_()
  //initializationReactor_(NULL)
{
  //PhaseRank initPhase = parent->schedule()->correction(id.rank());
  //String activityName = String("NonHomogenousInit_").append(toString(id));
  
  //Activity * activity = activityManagerInstance()->activityNew(activityName).ptr();
  //activity->phaseIs(initPhase);
  //activity->iterationIs(activity->currentIteration());
  
  //initializationReactor_ = new InitializationReactor(activity, this, notifier);

  //activity->statusIs(Activity::scheduled);
}

void
NonHomogeneousBasisCollectorImpl::PropagationReactor::onFinalState() {
  if (nonHomogeneousState().vectorSize() != 0) {
    parent()->finalStateIs(sliceId(), notifier()->finalState());
  } else {
    nonHomogeneousState_ = notifier()->finalState();
  }
}

/*typedef NonHomogeneousBasisCollectorImpl Nhbci;

Nhbci::PropagationReactor::InitializationReactor::InitializationReactor(
    Activity * notifier,
    PropagationReactor * parent,
    DynamPropagator * nonHomogeneousPropagator) :
  Activity::Notifiee(notifier),
  parent_(parent),
  nonHomogeneousPropagator_(nonHomogeneousPropagator)
{}

void
Nhbci::PropagationReactor::InitializationReactor::onStatus() {
  if (notifier()->status() == Activity::executing) {
    size_t vectorSize = nonHomogeneousPropagator_->vectorSize();
    nonHomogeneousPropagator_->initialStateIs(DynamState(vectorSize, 0.0)); 
    parent_->nonHomogeneousState_ = nonHomogeneousPropagator_->finalState();
    parent_->notifierIs(nonHomogeneousPropagator_.ptr()); // Now allow notifications
  }
}*/

} /* end namespace Pita */ } /* end namespace Hts */
