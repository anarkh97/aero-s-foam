#include "LinearFineIntegratorManager.h"

namespace Pita {

LinearFineIntegratorManager::LinearFineIntegratorManager(LinearDynamOps::Manager * dom, const GeneralizedAlphaParameter & fp) :
  FineIntegratorManager(fp.timeStepSize()),
  dynamOpsManager_(dom),
  forwardParameter_(fp),
  backwardParameter_(GeneralizedAlphaParameter(-fp.timeStepSize(), fp.rhoInfinity()))
{}

const GeneralizedAlphaParameter &
LinearFineIntegratorManager::parameter(HalfTimeSlice::Direction direction) const {
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

LinearGenAlphaIntegrator *
LinearFineIntegratorManager::createFineIntegrator(HalfTimeSlice::Direction direction) const {
  switch (direction) {
    case HalfTimeSlice::NO_DIRECTION:
      return NULL;
      break;
    case HalfTimeSlice::FORWARD: // Fall through
    case HalfTimeSlice::BACKWARD:
      return new LinearGenAlphaIntegrator(dynamOpsManager(), parameter(direction));
      break;
  }

  throw Fwk::InternalException("In LinearFineIntegratorManager::createFineIntegrator");
}

} // end namespace Pita
