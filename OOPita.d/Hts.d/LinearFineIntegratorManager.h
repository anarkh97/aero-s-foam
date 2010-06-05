#ifndef PITA_HTS_LINEARFINEINTEGRATORMANAGER_H
#define PITA_HTS_LINEARFINEINTEGRATORMANAGER_H

#include "FineIntegratorManager.h"

#include "../LinearDynamOps.h"

namespace Pita { namespace Hts {

template <typename GenAlphaIntegratorType>
class LinearFineIntegratorManager : public GenFineIntegratorManager<GenAlphaIntegratorType> {
public:
  EXPORT_PTRINTERFACE_TYPES(LinearFineIntegratorManager);

  LinearDynamOps::Manager * dynamOpsManager() const { return dynamOpsManager_.ptr(); }
  const GeneralizedAlphaParameter & parameter(HalfTimeSlice::Direction direction) const;

  static Ptr New(LinearDynamOps::Manager * dynOpsMgr, const GeneralizedAlphaParameter & forwardParam) {
    return new LinearFineIntegratorManager(dynOpsMgr, forwardParam);
  }

protected:
  LinearFineIntegratorManager(LinearDynamOps::Manager * dom, const GeneralizedAlphaParameter & fp);

  virtual GenAlphaIntegratorType * createFineIntegrator(HalfTimeSlice::Direction direction) const; // Overriden

private:
  LinearDynamOps::Manager::Ptr dynamOpsManager_;
  GeneralizedAlphaParameter forwardParameter_;
  GeneralizedAlphaParameter backwardParameter_;
};

template <typename GenAlphaIntegratorType>
LinearFineIntegratorManager<GenAlphaIntegratorType>::LinearFineIntegratorManager(LinearDynamOps::Manager * dom, const GeneralizedAlphaParameter & fp) :
  GenFineIntegratorManager<GenAlphaIntegratorType>(fp.timeStepSize()),
  dynamOpsManager_(dom),
  forwardParameter_(fp),
  backwardParameter_(GeneralizedAlphaParameter(-fp.timeStepSize(), fp.rhoInfinity()))
{}

template <typename GenAlphaIntegratorType>
const GeneralizedAlphaParameter &
LinearFineIntegratorManager<GenAlphaIntegratorType>::parameter(HalfTimeSlice::Direction direction) const {
  switch (direction) {
    case HalfTimeSlice::NO_DIRECTION: // Fall through
    case HalfTimeSlice::FORWARD:
      return forwardParameter_;
      break;
    case HalfTimeSlice::BACKWARD:
      return backwardParameter_; 
      break;
  }

  throw Fwk::InternalException("In LinearFineIntegratorManager::parameter");
}

template <typename GenAlphaIntegratorType>
GenAlphaIntegratorType *
LinearFineIntegratorManager<GenAlphaIntegratorType>::createFineIntegrator(HalfTimeSlice::Direction direction) const {
  switch (direction) {
    case HalfTimeSlice::NO_DIRECTION:
      return NULL;
      break;
    case HalfTimeSlice::FORWARD: // Fall through
    case HalfTimeSlice::BACKWARD:
      return new GenAlphaIntegratorType(this->dynamOpsManager(), parameter(direction));
      break;
  }

  throw Fwk::InternalException("In LinearFineIntegratorManager::createFineIntegrator");
}


} /* end namespace Hts */ } /* end namespace Pita */

#endif /* PITA_HTS_HTS_LINEARFINEINTEGRATORMANAGER_H */
