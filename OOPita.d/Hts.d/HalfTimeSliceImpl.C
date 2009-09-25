#include "HalfTimeSliceImpl.h"

#include "../IntegratorPropagator.h"

#include <memory>

namespace Pita { namespace Hts {

/* HalfTimeSliceImpl */

HalfTimeSliceImpl::HalfTimeSliceImpl(HalfSliceRank r,
                                     HalfTimeSlice::Direction d,
                                     DynamPropagator * dp) :
  HalfTimeSlice(r, d),
  propagator_(dp),
  previousSeedState_(),
  schedulingReactor_(NULL),
  localPropagationReactor_(NULL)
{}

/*void
HalfTimeSliceImpl::phaseIs(PhaseRank p) {
  localPropagationReactor_->notifier()->phaseIs(p);
  setPhase(p);
}*/

void
HalfTimeSliceImpl::seedIs(const Seed * s) {
  previousSeedState_ = s->state();
  //schedulingReactor_->notifierIs(s);
  setSeed(s);
}

void
HalfTimeSliceImpl::iterationIs(IterationRank i) {
  propagateSeed();
  propagatedSeed()->iterationIs(i);
}

void
HalfTimeSliceImpl::propagatedSeedIs(Seed * ps) {
  // TODO: Update ps immediately ?
  setPropagatedSeed(ps);
}

void
HalfTimeSliceImpl::propagateSeed() {
  if (seed()) {
    if (previousSeedState_.vectorSize() != 0) {
      //log() << "Reuse state\n";
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

/* HalfTimeSlice helper class definitions */ 

HalfTimeSliceImpl::SchedulingReactor::SchedulingReactor(const Seed * n, Activity * a) :
  Seed::NotifieeConst(n),
  activity_(a)
{}

void
HalfTimeSliceImpl::SchedulingReactor::notifierIs(const Seed * s) {
  if (s != notifier()) {
    Seed::NotifieeConst::notifierIs(s);
    if (this->notifier() != NULL)
      this->onState();
  }
}

inline
void
HalfTimeSliceImpl::SchedulingReactor::scheduleIfCondition() {
  if (this->schedulingCondition() && activity()->status() != Activity::scheduled) {
    activity()->iterationIs(activity()->currentIteration());
    activity()->statusIs(Activity::scheduled);
  }
}

void
HalfTimeSliceImpl::SchedulingReactor::onState() {
  scheduleIfCondition();
}

void
HalfTimeSliceImpl::SchedulingReactor::onStatus() {
  scheduleIfCondition();
  // TODO Improve protocol when ACTIVE <=> CONVERGED ?
}

bool
HalfTimeSliceImpl::ForwardSchedulingReactor::schedulingCondition() const { 
  return notifier()->status() != Seed::INACTIVE;
}

bool
HalfTimeSliceImpl::BackwardSchedulingReactor::schedulingCondition() const {
  switch (notifier()->status()) {
    case Seed::INACTIVE:
      return false;
    case Seed::ACTIVE:
      return true;
    case Seed::CONVERGED:
       return false;
    case Seed::SPECIAL:
       return true;
     default:
       break;
  }
  return false;
}

HalfTimeSliceImpl::LocalPropagationReactor::LocalPropagationReactor(Activity * notifier, HalfTimeSliceImpl * hs) :
  Activity::Notifiee(notifier),
  halfSlice_(hs)
{}

void
HalfTimeSliceImpl::LocalPropagationReactor::onStatus() {
  if (halfSlice_) {
    switch (notifier()->status()) {
      case Activity::executing:
        halfSlice_->propagateSeed();
        break;
      default:
        break;
    }
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
  
  // Setup internals
  //String activityName = String("LocalPropagation_") + toString(HalfSliceId(id.rank(), id.direction()));
  //Activity::Ptr activity = activityManagerInstance()->activityNew(activityName);
  //activity->phaseIs(halfSlice->phase());

  /*SchedulingReactor::Ptr schedulingReactor;
  if (id.direction() == HalfTimeSlice::FORWARD) {
    schedulingReactor = new ForwardSchedulingReactor(halfSlice->seed(), activity.ptr());
  } else {
    schedulingReactor = new BackwardSchedulingReactor(halfSlice->seed(), activity.ptr());
  }
  LocalPropagationReactor::Ptr localPropagationReactor = new LocalPropagationReactor(activity.ptr(), halfSlice.get());*/

  //halfSlice->schedulingReactor_ = schedulingReactor;
  //halfSlice->localPropagationReactor_ = localPropagationReactor;

  return halfSlice.release();
}

} /* end namespace Hts */ } /* end namespace Pita */
