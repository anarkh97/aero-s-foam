#ifndef PITA_HTS_HALFTIMESLICEIMPL_H
#define PITA_HTS_HALFTIMESLICEIMPL_H

#include "HalfTimeSlice.h"
#include "../AffineDynamPropagator.h"

#include "Fwk.h"

namespace Pita { namespace Hts {

class AffineHalfTimeSliceImpl : public HalfTimeSlice {
public:
  EXPORT_PTRINTERFACE_TYPES(AffineHalfTimeSliceImpl);

  class Manager;

  // Overriden mutators
  virtual void seedIs(const Seed * s);
  virtual void propagatedSeedIs(Seed * ps);

  // Added attributes 
  const AffineDynamPropagator * propagator() const { return propagator_.ptr(); }
  AffineDynamPropagator * propagator() { return propagator_.ptr(); }

  // Overriden
  virtual void iterationIs(IterationRank i);

protected:
  AffineHalfTimeSliceImpl(HalfSliceRank r, HalfTimeSlice::Direction d, AffineDynamPropagator * propagator);

  void propagateSeed(); 

private:
  AffineDynamPropagator::Ptr propagator_;
  DynamState previousSeedState_;

  friend class Manager;
};

// Helper classes

class AffineHalfTimeSliceImpl::Manager : public HalfTimeSlice::Manager, private Fwk::GenManagerImpl<AffineHalfTimeSliceImpl, HalfSliceId> {
public:
  EXPORT_PTRINTERFACE_TYPES(AffineHalfTimeSliceImpl::Manager);

  typedef Fwk::GenManagerInterface<AffineDynamPropagator*, HalfSliceId> PropagatorManager;

  virtual AffineHalfTimeSliceImpl * instance(const HalfSliceId & id) const;
  virtual size_t instanceCount() const;

  virtual AffineHalfTimeSliceImpl * instanceNew(const HalfSliceId & id);
  virtual void instanceDel(const HalfSliceId & id);

  static Ptr New(PropagatorManager * propagatorManager) {
    return new Manager(propagatorManager);
  }

protected:
  explicit Manager(PropagatorManager * propagatorManager);

  virtual AffineHalfTimeSliceImpl * createNewInstance(const HalfSliceId & id);

private:
  typedef Fwk::GenManagerImpl<AffineHalfTimeSliceImpl, HalfSliceId> Impl;

  PropagatorManager::Ptr propagatorManager_;
};

} /* end namespace Hts */ } /* end namespace Pita */

#endif /* PITA_HTS_HALFTIMESLICEIMPL_H */
