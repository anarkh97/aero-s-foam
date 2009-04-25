#include "PostProcessingManager.h"

#include <Driver.d/GeoSource.h>
extern GeoSource * geoSource;

namespace Pita {

// PostProcessingReactor implementation

PostProcessingReactor::PostProcessingReactor(const LinearGenAlphaIntegrator * notifier, 
                                             LinearPostProcessor * target,
                                             PostProcessor::FileSetId outputFileSet) :
  LinearGenAlphaIntegrator::NotifieeConst(NULL),
  notifier_(NULL),
  target_(target),
  outputFileSet_(outputFileSet)
{
  this->notifierIs(notifier);
  if (!this->target() || this->target()->fileStatus(this->outputFileSet()) == PostProcessor::NO_FILE) {
    throw Fwk::RangeException("In PostProcessingReactor::PostProcessingReactor");
  }
}

void
PostProcessingReactor::notifierIs(const DynamTimeIntegrator * notifier) {
  // If notifier is not actually a LinearGenAlphaIntegrator, do not accept notifications
  this->notifier_ = dynamic_cast<const LinearGenAlphaIntegrator *>(notifier);
  LinearGenAlphaIntegrator::NotifieeConst::notifierIs(notifier_);
}

inline
void
PostProcessingReactor::performOutput() const {
  this->target_->outputNew(this->outputFileSet(), this->notifier_);
}

void
PostProcessingReactor::onInitialCondition() {
  // Reset file(s) before writing initial state
  this->target_->fileStatusIs(this->outputFileSet(), PostProcessor::CLOSED);
  this->target_->fileStatusIs(this->outputFileSet(), PostProcessor::OPEN);
  this->performOutput();
}

void
PostProcessingReactor::onCurrentCondition() {
  this->performOutput();
}

// PostProcessingManager implementation

PostProcessingManager::PostProcessingManager(SDDynamPostProcessor * basePostProcessor,
                                             int localFileCount,
                                             const int * localFileId) :
  postProcessor_(LinearPostProcessor::New(geoSource, localFileCount, localFileId, basePostProcessor))
{}

PostProcessor::FileSetId
PostProcessingManager::outputFileSet(const IntegratorPropagator * op) const {
  PropagatorReactorContainer::const_iterator it = propagatorReactor_.find(op);
  return it != propagatorReactor_.end() ? it->second->fileSet() : PostProcessor::FileSetId();
}

void
PostProcessingManager::outputFileSetIs(const IntegratorPropagator * op, PostProcessor::FileSetId fs) {
  PropagatorReactorContainer::iterator it = propagatorReactor_.lower_bound(op);
  if (it != propagatorReactor_.end() && it->first == op) {
    // The propagator is already associated with a FileSet
    if (fs != PostProcessor::FileSetId()) {
      // Update FileSet association
      it->second->fileSetIs(fs);
    } else {
      // Default value corresponds to no fileset, stop observing propagator 
      propagatorReactor_.erase(it);
    }
  } else {
    if (fs != PostProcessor::FileSetId()) {
      // Begin observing propagator with associated FileSet
      PropagatorReactor::Ptr reactor = new PropagatorReactor(op, this->postProcessor_.ptr(), fs);
      propagatorReactor_.insert(it, make_pair(op, reactor));
    }
  }
}

// PostProcessingManager::PropagatorReactor implementation

PostProcessingManager::PropagatorReactor::PropagatorReactor(const IntegratorPropagator * notifier,
                                                            LinearPostProcessor * target,
                                                            PostProcessor::FileSetId fileSet) :
  IntegratorPropagator::Notifiee(NULL),
  notifier_(notifier),
  target_(target),
  fileSet_(fileSet),
  integratorReactor_(NULL)
{
  this->notifierIs(notifier);
}

void
PostProcessingManager::PropagatorReactor::onInitialState() {
  const LinearGenAlphaIntegrator * integrator = dynamic_cast<const LinearGenAlphaIntegrator *>(this->notifier_->integrator());
  if (integrator) {
    // Setup the notifiee to collect the output data
    this->integratorReactor_ = PostProcessingReactor::New(integrator, this->target_, this->fileSet());
  }
}

void
PostProcessingManager::PropagatorReactor::onFinalState() {
  this->integratorReactor_ = NULL;
  this->target_->fileStatusIs(this->fileSet(), PostProcessor::CLOSED);
}

} // end namespace Pita
