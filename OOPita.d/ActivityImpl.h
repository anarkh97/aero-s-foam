#ifndef PITA_ACTIVITYIMPL_H
#define PITA_ACTIVITYIMPL_H

#include "Activity.h"

#include <map>
#include <queue>

#include "Activity.h"

namespace Pita {

Activity::Manager::Ptr activityManagerInstance();

namespace ActivityImpl {
    
/* Comparison class for activities */
/* Used as a LessThanComparator */
class ActivityComp : public std::binary_function<Activity::Ptr, Activity::Ptr, bool> {
public:
  ActivityComp() {}

  bool operator()(const Activity::Ptr & a, const Activity::Ptr & b) const {
    return (a->iteration() == b->iteration()) ? (a->phase() > b->phase()) : (a->iteration() > b->iteration());
  }
};
  
class ManagerImpl;

class ActivityImpl : public Activity {
protected:
  ActivityImpl(const String & name, Fwk::Ptr<ManagerImpl> manager, IterationRank iteration);
  
  virtual void statusIs(Status s);
  virtual void iterationIs(IterationRank i);
  virtual void phaseIs(PhaseRank p);

  virtual IterationRank currentIteration() const;
  virtual PhaseRank currentPhase() const;
  
  Fwk::Ptr<ManagerImpl> manager() const { return manager_; }

private:
  friend class ManagerImpl;
  Fwk::Ptr<ManagerImpl> manager_;
};

    
class ManagerImpl : public Activity::Manager {
public:
  typedef Fwk::Ptr<ManagerImpl> Ptr;
  typedef Fwk::Ptr<const ManagerImpl> PtrConst;

  virtual Activity::Ptr activityNew(const String & name);
  virtual Activity::Ptr activity(const String & name) const;
  virtual void activityDel(const String & name);
  virtual size_t activityCount() const { return activities_.size(); } 
  virtual size_t scheduledActivityCount() const { return scheduledActivities_.size(); } 
  
  virtual void lastActivityIs(const Activity::Ptr & activity);

  virtual void targetPhaseIs(IterationRank iteration, PhaseRank phase);
  
  static Fwk::Ptr<Activity::Manager> Instance();

private:  
  // Data members
  std::priority_queue<Activity::Ptr, std::vector<Activity::Ptr>, ActivityComp> scheduledActivities_;
  std::map<String, Activity::Ptr> activities_; // Pool of all activities

  // Singleton instance
  static Activity::Manager::Ptr activityInstance_;	
};
    
} // end namespace ActivityImpl

} // end namespace Pita

#endif /* PITA_ACTIVITYIMPL_H */

