#ifndef PITA_HTS_NLLOCALNETWORK_H
#define PITA_HTS_NLLOCALNETWORK_H

#include "Fwk.h"
#include "Types.h"

#include "LocalNetwork.h"

#include "SliceMapping.h"

#include "../Seed.h"
#include "../RemoteStateMpiImpl.h"

#include "HalfTimeSlice.h"
#include "../JumpBuilder.h"
#include "CorrectionReductor.h"
#include "CorrectionReconstructor.h"
#include "../SeedUpdater.h"

#include "NlProjectionNetwork.h"

#include "JumpConvergenceEvaluator.h"
#include "../SeedDifferenceEvaluator.h"

#include "LocalNetworkImpl.h"

#include <map>
#include <deque>

namespace Pita { namespace Hts {

using namespace LocalNetworkImpl;

class NlLocalNetwork : public LocalNetwork {
public:
  EXPORT_PTRINTERFACE_TYPES(NlLocalNetwork);

  virtual void statusIs(Status s);
  void applyConvergenceStatus(); // TODO: bad naming

  NlLocalNetwork(SliceMapping * mapping,
                 RemoteState::MpiManager * commMgr,
                 HalfTimeSlice::Manager * htsMgr,
                 CorrectionReductor::Manager * corrRedMgr,
                 CorrectionReconstructor::Manager * corrReconMgr,
                 BasisCondensationManager * condensMgr,
                 ProjectionBuildingFactory * projBuildMgr,
                 JumpConvergenceEvaluator * jumpCvgMgr,
                 NonLinSeedDifferenceEvaluator::Manager * jumpEvalMgr);

  MainSeedMap activeSeeds() const { return seeds_[activeParity()]; }
  
  TaskList activeFinePropagators() const { return mapToDeque(finePropagators_[activeParity()]); }
  TaskList activePropagatedSeedSyncs() const { return mapToDeque(propagatedSeedSyncs_[activeParity()]); }
  TaskList activeJumpBuilders() const { return mapToDeque(jumpBuilders_[activeParity()]); }
  TaskList activeCorrectionPropagators() const { return mapToDeque(correctionPropagators_[activeParity()]); }
  TaskList activeSeedUpdaters() const { return mapToDeque(seedUpdaters_[activeParity()]); }

  TaskList activeCondensations() const { return mapToDeque(condensations_[activeParity()]); }
  TaskList activeProjectionBuilders() const { return mapToDeque(projectionBuilders_[activeParity()]); }

protected:
  HalfTimeSlice::Manager * htsMgr() { return htsMgr_.ptr(); }
  JumpBuilder::Manager * jumpMgr() { return jumpMgr_.ptr(); }
  NonLinSeedDifferenceEvaluator::Manager * jumpEvalMgr() { return jumpEvalMgr_.ptr(); }
  SeedUpdater::Manager * seedUpMgr() { return seedUpMgr_.ptr(); }
  CorrectionReductor::Manager * corrRedMgr() { return corrRedMgr_.ptr(); }
  CorrectionReconstructor::Manager * corrReconMgr() { return corrReconMgr_.ptr(); }
  BasisCondensationManager * condensMgr() { return condensMgr_.ptr(); }
  ProjectionBuildingFactory * projBuildMgr() { return projBuildMgr_.ptr(); }

  void init();

  void addMainSeed(HalfSliceRank seedRank);

  void addForwardPropagation(HalfSliceRank mainSeedRank);
  void addBackwardPropagation(HalfSliceRank mainSeedRank);

  void addForwardCondensation(HalfSliceRank mainSeedRank);
  void addBackwardCondensation(HalfSliceRank mainSeedRank);

  void addPropagatedSeedSend(HalfSliceRank seedRank);
  void addPropagatedSeedRecv(HalfSliceRank seedRank);

  void addJumpBuilder(HalfSliceRank seedRank);
 
  void addProjectionBuilding(HalfSliceRank seedRank);

  void addCorrectionReductor(HalfSliceRank seedRank);
  void addCorrectionReconstructor(HalfSliceRank seedRank);
  void addCorrectionSend(HalfSliceRank seedRank);
  void addCorrectionRecv(HalfSliceRank seedRank);
  void addReducedCorrectionSend(HalfSliceRank seedRank);
  void addReducedCorrectionRecv(HalfSliceRank seedRank);

  void addSeedUpdater(HalfSliceRank seedRank);

  void addNoCorrection(HalfSliceRank seedRank);

  NamedTask::Ptr forwardHalfSliceNew(HalfSliceRank sliceRank);
  NamedTask::Ptr backwardHalfSliceNew(HalfSliceRank sliceRank);
  
  NamedTask::Ptr propagatedSeedSendNew(HalfSliceRank seedRank);
  NamedTask::Ptr propagatedSeedRecvNew(HalfSliceRank seedRank);
  NamedTask::Ptr jumpBuilderNew(HalfSliceRank seedRank);

  NamedTask::Ptr correctionReductorNew(HalfSliceRank seedRank);
  NamedTask::Ptr correctionReconstructorNew(HalfSliceRank seedRank);
  NamedTask::Ptr correctionSendNew(HalfSliceRank seedRank);
  NamedTask::Ptr correctionRecvNew(HalfSliceRank seedRank);
  NamedTask::Ptr reducedCorrectionSendNew(HalfSliceRank seedRank);
  NamedTask::Ptr reducedCorrectionRecvNew(HalfSliceRank seedRank);

  NamedTask::Ptr seedUpdater(HalfSliceRank seedRank);
  NamedTask::Ptr seedUpdaterNew(HalfSliceRank seedRank);

private:
  HalfTimeSlice::Manager::Ptr htsMgr_;
  JumpBuilder::Manager::Ptr jumpMgr_;
  SeedUpdater::Manager::Ptr seedUpMgr_;
  CorrectionReductor::Manager::Ptr corrRedMgr_;
  CorrectionReconstructor::Manager::Ptr corrReconMgr_;

  BasisCondensationManager::Ptr condensMgr_;
  ProjectionBuildingFactory::Ptr projBuildMgr_;

  JumpConvergenceEvaluator::Ptr jumpCvgMgr_;

  NonLinSeedDifferenceEvaluator::Manager::Ptr jumpEvalMgr_;
  NoCorrectionManager::Ptr noCorrectionMgr_;
  
  MainSeedMap seeds_[2];

  typedef std::map<ActivationRange, NamedTask::Ptr, ActivationRange::Comparator> TaskMap;
  TaskMap finePropagators_[2];
  TaskMap propagatedSeedSyncs_[2];
  TaskMap jumpBuilders_[2];
  TaskMap correctionPropagators_[2];
  TaskMap seedUpdaters_[2];
  TaskMap condensations_[2];
  TaskMap projectionBuilders_[2];

  DISALLOW_COPY_AND_ASSIGN(NlLocalNetwork);
};

} /* end namespace Hts */ } /* end namespace Pita */

#endif /* PITA_HTS_NLLOCALNETWORK_H */
