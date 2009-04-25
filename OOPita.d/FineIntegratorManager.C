#include "FineIntegratorManager.h"

namespace Pita {

FineIntegratorManager::FineIntegratorManager(Seconds fineTimeStepSize) :
  fineTimeStepSize_(fineTimeStepSize),
  forwardFineIntegrator_(NULL),
  backwardFineIntegrator_(NULL)
{}

FineIntegratorManager::IntegratorType *
FineIntegratorManager::fineIntegrator(HalfTimeSlice::Direction direction) const {
  IntegratorType * integrator;
  
  switch (direction) {
    case HalfTimeSlice::FORWARD:
      if (!forwardFineIntegrator_) {
        forwardFineIntegrator_ = createFineIntegrator(direction);
      }
      integrator = forwardFineIntegrator_.ptr();
      break;

    case HalfTimeSlice::BACKWARD:
      if (!backwardFineIntegrator_) {
        backwardFineIntegrator_ = createFineIntegrator(direction);
      }
      integrator = backwardFineIntegrator_.ptr();
      break;

    default:
      integrator = NULL;
  }

  return integrator;
}

} // end namespace Pita
