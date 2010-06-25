#ifndef PITA_HTS_NONHOMOGENEOUSBASISCOLLECTORIMPL_H
#define PITA_HTS_NONHOMOGENEOUSBASISCOLLECTORIMPL_H

#include "BasisCollectorImpl.h"

#include "../SeedInitializer.h"

namespace Pita { namespace Hts {

class NonHomogeneousBasisCollectorImpl : public BasisCollectorImpl {
public:
  EXPORT_PTRINTERFACE_TYPES(NonHomogeneousBasisCollectorImpl);

  // New members
  DynamState nonHomogeneousState(const HalfSliceId & sliceId) const;

  static Ptr New() {
    return new NonHomogeneousBasisCollectorImpl();
  }

protected:
  class PropagationReactor;

  virtual BasisCollectorImpl::PropagationReactor * propagationReactorNew(
      const DynamPropagator * notifier,
      const HalfSliceId & id); // Overriden

  NonHomogeneousBasisCollectorImpl();
};

class NonHomogeneousBasisCollectorImpl::PropagationReactor :
  public BasisCollectorImpl::PropagationReactor {
public:
  EXPORT_PTRINTERFACE_TYPES(PropagationReactor);
  
  // Overriden members
  virtual void onFinalState();

  // Added members
  DynamState nonHomogeneousState() const { return nonHomogeneousState_; }

  PropagationReactor(const DynamPropagator * notifier,
                     const HalfSliceId & id,
                     NonHomogeneousBasisCollectorImpl * parent);

private:
  DynamState nonHomogeneousState_;
};

} /* end namespace Pita */ } /* end namespace Hts */

#endif /* PITA_HTS_NONHOMOGENEOUSBASISCOLLECTORIMPL_H */
