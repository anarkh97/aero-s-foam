#ifndef PITA_ACTIVITY_H
#define PITA_ACTIVITY_H

#include "Fwk.h"
#include "Types.h"

namespace Pita {

/* Activity sub-module to share the processor time between the differents objects of the Pita module
** Especially to handle dependencies between the different phases of the Pita algorithm
** The manager runs the activities until the scheduled queue becomes empty */
  
/* Auxiliary classes to determine priority */
class IterationRank;
class PhaseRank;

/* Activity interface class */

class Activity : public Fwk::PtrInterface<Activity> {
public:
  typedef Fwk::Ptr<Activity> Ptr;
  typedef Fwk::Ptr<const Activity> PtrConst;
    
  /* Notifiee class for Activities */
  class Notifiee : public Fwk::BaseNotifiee<Activity> {
  public:
	  typedef Fwk::Ptr<Notifiee> Ptr;
	  typedef Fwk::Ptr<const Notifiee> PtrConst;

    explicit Notifiee(Activity * act) : Fwk::BaseNotifiee<Activity>(act) {}

    virtual void onStatus() {}
    virtual void onIteration() {}
    virtual void onPhase() {}
  };
  
  Notifiee::Ptr lastNotifiee() const { return notifiee_; }
  void lastNotifieeIs(Notifiee * n) { notifiee_ = n; }

  class Manager;

  enum Status {
      free = 0,
      scheduled,
      executing
  };

  String name() const { return name_; }
  Status status() const { return status_; }
  IterationRank iteration() const { return iteration_; }
  virtual IterationRank currentIteration() const = 0;
  PhaseRank phase() const { return phase_; }
  virtual PhaseRank currentPhase() const = 0;
  
  virtual void statusIs(Status s) = 0;
  virtual void iterationIs(IterationRank rank) = 0;
  virtual void phaseIs(PhaseRank rank) = 0;

protected:
  explicit Activity(const String & name, IterationRank iteration, PhaseRank phase) :
    name_(name),
    notifiee_(NULL),
    status_(free),
    iteration_(iteration),
    phase_()
  {}

  void setStatus(Status s) {
    status_ = s;
    if (lastNotifiee())
      lastNotifiee()->onStatus();
  }

  void setIteration(IterationRank i) {
    iteration_ = i;
    if (lastNotifiee())
      lastNotifiee()->onIteration();
  }

  void setPhase(PhaseRank p) {
    phase_ = p;
    if (lastNotifiee())
      lastNotifiee()->onPhase();
  }

private:
  String name_;
  Notifiee * notifiee_;
  Status status_;
  IterationRank iteration_;
  PhaseRank phase_;
};


/* Manager interface */

class Activity::Manager : public Fwk::PtrInterface<Activity::Manager> {
public:
  typedef Fwk::Ptr<Activity::Manager> Ptr;

  virtual Fwk::Ptr<Activity> activityNew(const String & name) = 0;
  virtual Fwk::Ptr<Activity> activity(const String & name) const = 0;
  virtual void activityDel(const String & name) = 0;
  virtual size_t activityCount() const = 0;
  virtual size_t scheduledActivityCount() const = 0; 
  
  virtual void lastActivityIs(const Activity::Ptr & activity) = 0;

  virtual void targetPhaseIs(IterationRank iteration, PhaseRank phase) = 0; 
  
  IterationRank currentIteration() const { return currentIteration_; }  
  PhaseRank currentPhase() const { return currentPhase_; }

protected:
  Manager() : currentIteration_(), currentPhase_() {}
  
  void setCurrentIteration(IterationRank i) {
    currentIteration_ = i;
  }

  void setCurrentPhase(PhaseRank p) {
    currentPhase_ = p;
  }

private:
  IterationRank currentIteration_;
  PhaseRank currentPhase_;
};

/* Entry point */

extern Activity::Manager::Ptr activityManagerInstance();

} // end namespace Pita

#endif
