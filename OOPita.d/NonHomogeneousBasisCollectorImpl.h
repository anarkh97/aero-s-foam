#ifndef PITA_HTS_NONHOMOGENEOUSBASISCOLLECTORIMPL_H
#define PITA_HTS_NONHOMOGENEOUSBASISCOLLECTORIMPL_H

#include "HalfSliceBasisCollectorImpl.h"

#include "SeedInitializer.h"
#include "HalfSliceSchedule.h"

#include "Activity.h"

namespace Pita { namespace Hts {

class NonHomogeneousBasisCollectorImpl : public HalfSliceBasisCollectorImpl {
public:
  EXPORT_PTRINTERFACE_TYPES(NonHomogeneousBasisCollectorImpl);

  // New members
  DynamState nonHomogeneousState(const HalfSliceId & sliceId) const;

  static Ptr New() {
    return new NonHomogeneousBasisCollectorImpl();
  }

protected:
  class PropagationReactor;

  virtual HalfSliceBasisCollectorImpl::PropagationReactor * propagationReactorNew(
      IntegratorPropagator * notifier,
      const HalfSliceId & id); // Overriden

  NonHomogeneousBasisCollectorImpl();
};

class NonHomogeneousBasisCollectorImpl::PropagationReactor :
  public HalfSliceBasisCollectorImpl::PropagationReactor {
public:
  EXPORT_PTRINTERFACE_TYPES(PropagationReactor);
  
  // Overriden members
  virtual void onFinalState();

  // Added members
  DynamState nonHomogeneousState() const { return nonHomogeneousState_; }

  PropagationReactor(IntegratorPropagator * notifier,
                     const HalfSliceId & id,
                     NonHomogeneousBasisCollectorImpl * parent);

protected:
  //class InitializationReactor;
  //friend class InitializationReactor;

private:
  DynamState nonHomogeneousState_;

  //Fwk::Ptr<InitializationReactor> initializationReactor_;
};

/*class NonHomogeneousBasisCollectorImpl::PropagationReactor::InitializationReactor :
  public Activity::Notifiee {
public:
  EXPORT_PTRINTERFACE_TYPES(InitializationReactor);

  virtual void onStatus(); // Overriden

  InitializationReactor(Activity * notifier,
                        PropagationReactor * parent,
                        DynamPropagator * nonHomogeneousPropagator);

private:
  NonHomogeneousBasisCollectorImpl::PropagationReactor * parent_;
  DynamPropagator::Ptr nonHomogeneousPropagator_;
};*/

} /* end namespace Pita */ } /* end namespace Hts */

#endif /* PITA_HTS_NONHOMOGENEOUSBASISCOLLECTORIMPL_H */
