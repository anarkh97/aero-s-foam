#ifndef PITA_HTS_NLTASKMANAGER_H
#define PITA_HTS_NLTASKMANAGER_H

#include "../TaskManager.h"

#include "SliceMapping.h"
#include "NlPropagatorManager.h"
#include "../SeedInitializer.h"
#include "../PostProcessingManager.h"
#include "../SeedErrorEvaluator.h"
#include "../RemoteStateMpiImpl.h"
#include "NlProjectionNetwork.h"
#include "PropagationDataSharing.h"

#include "NlLocalNetwork.h"

namespace Pita { namespace Hts {

class NlTaskManager : public TaskManager {
public:
  EXPORT_PTRINTERFACE_TYPES(NlTaskManager);

  virtual void iterationInc(); // overriden
  virtual Phase * phase();     // overriden
  virtual void phaseInc();     // overriden

  CpuRank localCpu() const { return commMgr_->localCpu(); }

  NlTaskManager(SliceMapping * mapping, RemoteState::MpiManager * commMgr,
                NlPropagatorManager * propagatorMgr, SeedInitializer * seedInitializer,
                PostProcessing::Manager * postProcessingMgr, SeedErrorEvaluator::Manager * jumpEvaluatorMgr,
                double projectorTolerance, IterationRank lastIteration);

protected:
  void initialize();

  // Execution
  void setPhase(Phase * p) { phase_ = p; }

  typedef void (NlTaskManager::*Continuation)();
  void setContinuation(Continuation c) { continuation_ = c; }

  void noop() {}
  void scheduleNothing();

  void scheduleInitialSeed();
  void fillFirstProjectionBasis(); 
  void scheduleProjectionBasisCondensation();
  void scheduleFinePropagation();
  void scheduleDataSharing();
  void checkConvergence();
  void schedulePropagatedSeedSynchronization();
  void scheduleJumpProjection();
  void scheduleProjectionBuilding();
  void scheduleCorrectionPropagation();
  void scheduleSeedUpdate();
  void enrichProjectionBasis();

private:
  SliceMapping::Ptr mapping_;
  RemoteState::MpiManager::Ptr commMgr_;

  NlPropagatorManager::Ptr propagatorMgr_;
  SeedInitializer::Ptr seedInitializer_;
  PostProcessing::Manager::Ptr postProcessingMgr_;
  SeedErrorEvaluator::Manager::Ptr jumpEvaluatorMgr_;

  double projectorTolerance_;

  Phase::Ptr phase_;
  Continuation continuation_;

  NlLocalNetwork::Ptr localNetwork_;
  
  PropagationDataSharing::Ptr sharing_;
  NlProjectionNetwork::Ptr projectionNetwork_;
  IterationRank lastIteration_;
};

} /* end namespace Hts */ } /* end namespace Pita */

#endif /* PITA_HTS_NLTASKMANAGER_H */
