#ifndef PITA_HTS_LINEARLOCALNETWORK_H
#define PITA_HTS_LINEARLOCALNETWORK_H

#include "Fwk.h"
#include "Types.h"

#include "SliceMapping.h"

#include "../Seed.h"
#include "../RemoteState.h"

#include "HalfTimeSlice.h"
#include "../JumpProjector.h"
#include "../UpdatedSeedAssembler.h"
#include "../CorrectionPropagator.h"

#include "../SeedDifferenceEvaluator.h"

#include <map>
#include <deque>
#include <cassert>

#include "LocalNetwork.h"
#include "LocalNetworkImpl.h"

namespace Pita { namespace Hts {

using namespace LocalNetworkImpl;

class ReducedCorrectionManager : public Fwk::PtrInterface<ReducedCorrectionManager> {
public:
  EXPORT_PTRINTERFACE_TYPES(ReducedCorrectionManager);

  JumpProjector::Manager * jumpProjMgr() const { return jumpProjMgr_.ptr(); }
  CorrectionPropagator<Vector>::Manager * rcpMgr() const { return rcpMgr_.ptr(); }
  CorrectionPropagator<DynamState>::Manager * fcpMgr() const { return fcpMgr_.ptr(); }
  UpdatedSeedAssembler::Manager * usaMgr() const { return usaMgr_.ptr(); }
  
  ReducedCorrectionManager(JumpProjector::Manager * jumpProjMgr,
                           CorrectionPropagator<Vector>::Manager * rcpMgr,
                           CorrectionPropagator<DynamState>::Manager * fcpMgr,  // can be NULL
                           UpdatedSeedAssembler::Manager * usaMgr) :
    jumpProjMgr_(jumpProjMgr),
    rcpMgr_(rcpMgr),
    fcpMgr_(fcpMgr),
    usaMgr_(usaMgr)
  {}

private:
  JumpProjector::Manager::Ptr jumpProjMgr_;
  CorrectionPropagator<Vector>::Manager::Ptr rcpMgr_;
  CorrectionPropagator<DynamState>::Manager::Ptr fcpMgr_;
  UpdatedSeedAssembler::Manager::Ptr usaMgr_;

  DISALLOW_COPY_AND_ASSIGN(ReducedCorrectionManager);
};


class LinearLocalNetwork : public LocalNetwork {
public:
  EXPORT_PTRINTERFACE_TYPES(LinearLocalNetwork);
 
  virtual void statusIs(Status s) { /* TODO */ }
  virtual void convergedSlicesInc(); // overriden

  /* Local network elements */
  TaskList halfTimeSlices() const;
  TaskList activeHalfTimeSlices() const;
  TaskList activeJumpProjectors() const;
  TaskList activeLeftSeedSyncs() const;
  TaskList activeSeedAssemblers() const;

  TaskList activeFullTimeSlices() const;
  TaskList activeCorrectionSyncs() const;

  TaskList activeCoarseTimeSlices() const;
  TaskList activeFullCorrectionSyncs() const;
  
  SeedMap mainSeeds() const;
  MainSeedMap activeMainSeeds() const;

  LinearLocalNetwork(SliceMapping * mapping,
                     HalfTimeSlice::Manager * sliceMgr,
                     ReducedCorrectionManager * redCorrMgr,
                     RemoteState::Manager * commMgr,
                     LinSeedDifferenceEvaluator::Manager * jumpErrorMgr);

protected:
  void init();

  void buildForwardPropagation(HalfSliceRank sliceRank);
  void buildBackwardPropagation(HalfSliceRank sliceRank);
  void buildPrimalCorrectionNetwork(HalfSliceRank sliceRank);
  void buildDualCorrectionNetwork(HalfSliceRank sliceRank);
  
  virtual void buildCorrectionPropagator(HalfSliceRank sliceRank);
  virtual void buildCorrectionSynchronizationSend(HalfSliceRank sliceRank);
  virtual void buildCorrectionSynchronizationRecv(HalfSliceRank sliceRank);
  
  void buildPropagatedSeedSend(HalfSliceRank sliceRank);
  void buildPropagatedSeedRecv(HalfSliceRank sliceRank);
  void buildJumpEstimator(HalfSliceRank sliceRank);
  void buildJumpBuilder(HalfSliceRank sliceRank);
  void buildSeedUpdater(HalfSliceRank sliceRank);
  
  void buildReducedCorrectionPropagator(HalfSliceRank sliceRank);
  void buildFullCorrectionPropagator(HalfSliceRank sliceRank);
  void buildReducedCorrectionSynchronization(HalfSliceRank sliceRank);
  void buildFullCorrectionSynchronization(HalfSliceRank sliceRank);

  void eraseInactive(TaskMap & task) { task.erase(task.begin(), task.lower_bound(firstActiveSlice())); }
  TaskList getActive(const TaskMap task[2]) const;
  static TaskList getAll(const TaskMap & task);

private:
  HalfTimeSlice::Manager::Ptr sliceMgr_;
  ReducedCorrectionManager::Ptr redCorrMgr_;
  
  TaskMap halfTimeSlice_[2];
  TaskMap jumpProjector_[2];
  TaskMap seedAssembler_[2];
  TaskMap leftSeedSync_[2];

  TaskMap correctionPropagator_[2];
  TaskMap correctionSync_[2];

  TaskMap alternateCorrectionPropagator_[2];
  TaskMap alternateCorrectionSync_[2];
  
  SeedMap mainSeed_; 

  LinSeedDifferenceEvaluator::Manager::Ptr jumpErrorMgr_;
  
  CorrectionPropagatorBuilder<DynamState> fullCorrectionBuilder_;
  CorrectionPropagatorBuilder<Vector> reducedCorrectionBuilder_;

  NoCorrectionManager::Ptr noCorrectionMgr_;

  DISALLOW_COPY_AND_ASSIGN(LinearLocalNetwork);
};

} /* end namespace Hts */ } /* end namespace Pita */

#endif /* PITA_HTS_LINEARLOCALNETWORK_H */
