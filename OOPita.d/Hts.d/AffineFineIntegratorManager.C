#include "AffineFineIntegratorManager.h"

namespace Pita { namespace Hts {

AffineFineIntegratorManager::AffineFineIntegratorManager(
    LinearDynamOps::Manager * dom,
    const GeneralizedAlphaParameter & fp,
    const Schedule * schedule) : 
  LinearFineIntegratorManager<AffineGenAlphaIntegrator>(dom, fp),
  schedule_(schedule),
  reactor_()
{}
  
AffineGenAlphaIntegrator *
AffineFineIntegratorManager::createFineIntegrator(HalfTimeSlice::Direction direction) const {
  AffineGenAlphaIntegrator * integrator = LinearFineIntegratorManager<AffineGenAlphaIntegrator>::createFineIntegrator(direction);
  if (integrator) {
    Activity::Ptr activity = activityManagerInstance()->activityNew(String("AffineForce_") + toString(direction));
    activity->iterationIs(IterationRank(0));
    activity->phaseIs(schedule()->endIteration());
    activity->statusIs(Activity::scheduled);

    reactor_.push_back(AffineForceReactor::New(activity.ptr(), integrator));
  }   
  return integrator;
}

AffineFineIntegratorManager::AffineForceReactor::AffineForceReactor(
    Activity * notifier,
    AffineGenAlphaIntegrator * target) : 
  Activity::Notifiee(notifier),
  target_(target)
{}

void
AffineFineIntegratorManager::AffineForceReactor::onStatus() {
  if (notifier()->status() == Activity::executing) {
    target_->externalForceFlagIs(false);
  }
}

} /* end namespace Hts */ } /* end namespace Pita */
