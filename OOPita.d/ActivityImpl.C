#include "ActivityImpl.h"
#include <iostream>

#include <Timers.d/GetTime.h>

namespace Pita {

Fwk::Ptr<Activity::Manager> activityManagerInstance() {
    return ActivityImpl::ManagerImpl::Instance();
}

namespace ActivityImpl {

// Definition of static member
Fwk::Ptr<Activity::Manager> ManagerImpl::activityInstance_ = NULL;

// ActivityImpl implementation

ActivityImpl::ActivityImpl(const String & name, Fwk::Ptr<ManagerImpl> manager, IterationRank iteration) :
  Activity(name, iteration, PhaseRank()),
  manager_(manager)
{}

void
ActivityImpl::statusIs(Status s) {
  if (s == status())
    return;

  if ((status() == free && s == executing) || (status() == scheduled && s == free)) {
    throw Fwk::RangeException("In ActivityImpl::statusIs"); // Unhandled transitions -- TODO Improve
  } else if (s == scheduled) {
    manager_->lastActivityIs(this);
  }

  setStatus(s);
}

void
ActivityImpl::iterationIs(IterationRank i) {
  if (iteration() != i)
    setIteration(i);
}

void
ActivityImpl::phaseIs(PhaseRank p) {
  if (phase() != p)
    setPhase(p);
}
  
IterationRank
ActivityImpl::currentIteration() const {
  return manager_->currentIteration();
}

PhaseRank
ActivityImpl::currentPhase() const {
  return manager_->currentPhase();
}


// ManagerImpl implementation

Activity::Manager::Ptr
ManagerImpl::Instance() {
  if (activityInstance_ == NULL) {
    activityInstance_ = new ManagerImpl();
  }
  return activityInstance_;
}
  
Activity::Ptr
ManagerImpl::activityNew(const String & name) {
  Activity::Ptr activity = activities_[name];
  if (activity != NULL) {
    throw Fwk::NameInUseException();
  }
  activity = new ActivityImpl(name, this, currentIteration());
  activities_[name] = activity;
  return activity;
}

Activity::Ptr
ManagerImpl::activity(const String & name) const {
  std::map<String, Activity::Ptr>::const_iterator it = activities_.find(name);
  if (it != activities_.end() ) {
    return (*it).second;
  }
  return NULL; 
}
  
void
ManagerImpl::activityDel(const String & name) {
  activities_.erase(name);
}
  
void
ManagerImpl::lastActivityIs(const Activity::Ptr & activity) {
  //log() << "(" << currentIteration() << ", " << currentPhase() << ") ";
  //log() << "Scheduling " << activity->name() << " at (" << activity->iteration() << ", " << activity->phase() << ")\n";
  scheduledActivities_.push(activity);
}

void
ManagerImpl::targetPhaseIs(IterationRank targetIteration, PhaseRank targetPhase) {
  // Find the most recent activites to run and run them in order
  while (!scheduledActivities_.empty()) {
    // Determine the next activity to run
    Activity::Ptr nextToRun = scheduledActivities_.top();
    // If the next time is greater than the specified time, break
    IterationRank nextIteration = nextToRun->iteration();
    PhaseRank nextPhase = nextToRun->phase();
    if (nextIteration > targetIteration || (nextIteration == targetIteration && nextPhase > targetPhase)) {
      break;
    }
    setCurrentIteration(nextIteration);
    setCurrentPhase(nextPhase);
    
    // Run the most urgent activity and remove it from the queue
    scheduledActivities_.pop();
    
    double tic = getTime();

    //log() << "(" << currentIteration() << ", " << currentPhase() << ") ";
    //log() << "Running " << nextToRun->name() << "\n";

    nextToRun->statusIs(Activity::executing);
    
    double toc = getTime();
    log() << nextToRun->name() << " took " << (toc - tic) / 1000.0 << " s\n";
    
    nextToRun->statusIs(Activity::free);
  }
  setCurrentIteration(targetIteration);
  setCurrentPhase(targetPhase);
}

} // end namespace ActivityImpl

} // end namepsace Pita
