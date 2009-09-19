#ifndef PITA_HALFTIMESLICEIMPL_H
#define PITA_HALFTIMESLICEIMPL_H

#include "HalfTimeSlice.h"
#include "DynamPropagator.h"

#include "DynamTimeIntegrator.h" 

#include "Fwk.h"

namespace Pita {

class HalfTimeSliceImpl : public HalfTimeSlice {
public:
  typedef Fwk::Ptr<HalfTimeSliceImpl> Ptr;
  typedef Fwk::Ptr<const HalfTimeSliceImpl> PtrConst;

  class Manager;
  class SchedulingReactor;
  class ForwardSchedulingReactor;
  class BackwardSchedulingReactor;
  class LocalPropagationReactor;

  // Overriden mutators
  //virtual void phaseIs(PhaseRank p);
  virtual void seedIs(const Seed * s);
  virtual void propagatedSeedIs(Seed * ps);

  // Added attributes 
  const DynamPropagator * propagator() const { return propagator_.ptr(); }
  DynamPropagator * propagator() { return propagator_.ptr(); }

  // Overriden
  virtual void iterationIs(IterationRank i);

protected:
  HalfTimeSliceImpl(HalfSliceRank r,
                    HalfTimeSlice::Direction d,
                    DynamPropagator * propagator);
  void propagateSeed(); 

private:

  DynamPropagator::Ptr propagator_;
  DynamState previousSeedState_;

  Fwk::Ptr<SchedulingReactor> schedulingReactor_;
  Fwk::Ptr<LocalPropagationReactor> localPropagationReactor_;

  friend class Manager;
  friend class LocalPropagationReactor;
};

// Helper classes

class HalfTimeSliceImpl::Manager : public HalfTimeSlice::Manager, private Fwk::GenManagerImpl<HalfTimeSliceImpl, HalfSliceId> {
public:
  typedef Fwk::Ptr<HalfTimeSliceImpl::Manager> Ptr;
  typedef Fwk::Ptr<const HalfTimeSliceImpl::Manager> PtrConst;

  typedef Fwk::GenManagerInterface<DynamPropagator*, HalfSliceId> PropagatorManager;

  virtual HalfTimeSliceImpl * instance(const HalfSliceId & id) const;
  virtual size_t instanceCount() const;

  virtual HalfTimeSliceImpl * instanceNew(const HalfSliceId & id);
  virtual void instanceDel(const HalfSliceId & id);

  static Ptr New(PropagatorManager * propagatorManager) {
    return new Manager(propagatorManager);
  }

protected:
  explicit Manager(PropagatorManager * propagatorManager);

  virtual HalfTimeSliceImpl * createNewInstance(const HalfSliceId & id);

private:
  typedef Fwk::GenManagerImpl<HalfTimeSliceImpl, HalfSliceId> Impl;

  PropagatorManager::Ptr propagatorManager_;
};

class HalfTimeSliceImpl::SchedulingReactor : public Seed::NotifieeConst {
public:
  typedef Fwk::Ptr<SchedulingReactor> Ptr;
  typedef Fwk::Ptr<const SchedulingReactor> PtrConst;

  virtual void notifierIs(const Seed * s); // overriden

  virtual void onState(); // overriden
  virtual void onStatus(); // overriden

  Activity * activity() const { return activity_.ptr(); }

protected:
  SchedulingReactor(const Seed * notifier, Activity * activity);

  virtual bool schedulingCondition() const = 0;

private:
  void scheduleIfCondition();

  Activity::Ptr activity_;
};
  
class HalfTimeSliceImpl::LocalPropagationReactor : public Activity::Notifiee {
public:
  typedef Fwk::Ptr<LocalPropagationReactor> Ptr;
  typedef Fwk::Ptr<const LocalPropagationReactor> PtrConst;

  virtual void onStatus(); // overriden

  LocalPropagationReactor(Activity * notifier, HalfTimeSliceImpl * hs);

private:
  HalfTimeSliceImpl * halfSlice_;
};

class HalfTimeSliceImpl::ForwardSchedulingReactor : public HalfTimeSliceImpl::SchedulingReactor {
public:
  ForwardSchedulingReactor(const Seed * notifier, Activity * activity) :
    HalfTimeSliceImpl::SchedulingReactor(notifier, activity)
  {}

protected:
  virtual bool schedulingCondition() const; // overriden
};

class HalfTimeSliceImpl::BackwardSchedulingReactor : public HalfTimeSliceImpl::SchedulingReactor {
public:
  BackwardSchedulingReactor(const Seed * notifier, Activity * activity) :
    HalfTimeSliceImpl::SchedulingReactor(notifier, activity)
  {}

protected:
  virtual bool schedulingCondition() const; // overriden
};

} // end namespace Pita

#endif /* PITA_HALFTIMESLICEIMPL_H */
