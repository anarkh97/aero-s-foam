#ifndef PITA_HTS_BASICHALFTIMESLICEIMPL_H
#define PITA_HTS_BASICHALFTIMESLICEIMPL_H

#include "HalfTimeSlice.h"

#include "../DynamPropagator.h"

#include "Fwk.h"

namespace Pita { namespace Hts {

class BasicHalfTimeSliceImpl : public HalfTimeSlice {
public:
  EXPORT_PTRINTERFACE_TYPES(BasicHalfTimeSliceImpl);

  class Manager;

  // Added attributes 
  const DynamPropagator * propagator() const { return propagator_.ptr(); }
  DynamPropagator * propagator() { return propagator_.ptr(); }

  // Overriden
  virtual void iterationIs(IterationRank i);

protected:
  BasicHalfTimeSliceImpl(HalfSliceRank r, HalfTimeSlice::Direction d, DynamPropagator * propagator);

  void propagateSeed();

private:
  DynamPropagator::Ptr propagator_;

  friend class Manager;
};

// Helper classes

class BasicHalfTimeSliceImpl::Manager : public HalfTimeSlice::Manager, private Fwk::GenManagerImpl<BasicHalfTimeSliceImpl, HalfSliceId> {
public:
  EXPORT_PTRINTERFACE_TYPES(Manager);

  typedef Fwk::GenManagerInterface<DynamPropagator*, HalfSliceId> DynamPropagatorManager;

  virtual BasicHalfTimeSliceImpl * instance(const HalfSliceId & id) const;
  virtual size_t instanceCount() const;

  virtual BasicHalfTimeSliceImpl * instanceNew(const HalfSliceId & id);
  virtual void instanceDel(const HalfSliceId & id);

  static Ptr New(DynamPropagatorManager * propagatorManager) {
    return new Manager(propagatorManager);
  }

protected:
  explicit Manager(DynamPropagatorManager * propagatorManager);

  virtual BasicHalfTimeSliceImpl * createNewInstance(const HalfSliceId & id);

private:
  typedef Fwk::GenManagerImpl<BasicHalfTimeSliceImpl, HalfSliceId> Impl;

  DynamPropagatorManager::Ptr propagatorManager_;
};

} /* end namespace Hts */ } /* end namespace Pita */

#endif /* PITA_HTS_BASICHALFTIMESLICEIMPL_H */
