#ifndef PITA_HTS_HALFTIMESLICEIMPL_H
#define PITA_HTS_HALFTIMESLICEIMPL_H

#include "HalfTimeSlice.h"
#include "../DynamPropagator.h"

#include "../DynamTimeIntegrator.h" 

#include "Fwk.h"

namespace Pita { namespace Hts {

class HalfTimeSliceImpl : public HalfTimeSlice {
public:
  EXPORT_PTRINTERFACE_TYPES(HalfTimeSliceImpl);

  class Manager;

  // Overriden mutators
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

  friend class Manager;
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

} /* end namespace Hts */ } /* end namespace Pita */

#endif /* PITA_HTS_HALFTIMESLICEIMPL_H */
