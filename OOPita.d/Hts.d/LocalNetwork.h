#ifndef PITA_HTS_LOCALNETWORK_H
#define PITA_HTS_LOCALNETWORK_H

#include "Fwk.h"
#include "Types.h"

#include "../TaskManager.h"
#include "../Seed.h"

#include "HalfTimeSlice.h"
#include "../JumpProjector.h"
#include "ReducedFullTimeSlice.h"
#include "../UpdatedSeedAssembler.h"

#include "../RemoteState.h"

#include <map>
#include <deque>

namespace Pita { namespace Hts {

class LocalNetwork : public Fwk::PtrInterface<LocalNetwork> {
public:
  EXPORT_PTRINTERFACE_TYPES(LocalNetwork);

  typedef std::map<HalfSliceRank, NamedTask::Ptr> TaskMap;
  typedef std::deque<NamedTask::Ptr> TaskList;

  typedef std::map<HalfSliceRank, Seed::Ptr> SeedMap;
  typedef std::map<SliceRank, Seed::Ptr> MainSeedMap;
  
  typedef std::deque<HalfTimeSlice::Ptr> HTSList;
  typedef std::deque<ReducedFullTimeSlice::Ptr> FTSList;
  
  /* Basic info */
  FullSliceCount totalFullSlices() const { return totalFullSlices_; }
  CpuCount availableCpus() const { return availableCpus_; }
  HalfSliceCount maxWorkload() const { return maxWorkload_; }
  CpuRank localCpu() const { return localCpu_; }

  /* Current state */
  HalfSliceRank firstActiveSlice() const { return HalfSliceRank(taskManager_->firstCurrentTask()); }
  HalfSliceRank firstInactiveSlice() const { return HalfSliceRank(taskManager_->firstWaitingTask()); }
  HalfSliceCount convergedSlices() const { return HalfSliceCount(taskManager_->completedTasks()); }
  void convergedSlicesInc();

  /* Local network elements */
  HTSList halfTimeSlices() const { return HTSList(halfTimeSlice_); }
  
  TaskList activeHalfTimeSlices() const;
  TaskList activeJumpProjectors() const;
  TaskList activeLeftSeedSyncs() const;
  TaskList activeFullTimeSlices() const;
  TaskList activeCorrectionSyncs() const;
  TaskList activeSeedAssemblers() const;

  MainSeedMap activeMainSeeds() const;

  LocalNetwork(FullSliceCount totalFullSlices,
               CpuCount availableCpus,
               HalfSliceCount maxWorkload,
               CpuRank localCpu,
               HalfTimeSlice::Manager * sliceMgr,
               JumpProjector::Manager * jumpProjMgr,
               ReducedFullTimeSlice::Manager * ftsMgr,
               UpdatedSeedAssembler::Manager * usaMgr,
               RemoteState::Manager * commMgr);

protected:
  void init();

  Seed * getSeed(const SeedId & id);
  ReducedSeed * getReducedSeed(const SeedId & id);

private:
  FullSliceCount totalFullSlices_;
  CpuCount availableCpus_;
  HalfSliceCount maxWorkload_;
  CpuRank localCpu_;

  TaskManager::Ptr taskManager_;
 
  HalfTimeSlice::Manager::Ptr sliceMgr_;
  JumpProjector::Manager::Ptr jumpProjMgr_;
  ReducedFullTimeSlice::Manager::Ptr ftsMgr_;
  UpdatedSeedAssembler::Manager::Ptr usaMgr_;

  Seed::Manager::Ptr seedMgr_;
  ReducedSeed::Manager::Ptr reducedSeedMgr_;

  RemoteState::Manager::Ptr commMgr_;
  
  HTSList halfTimeSlice_;
  TaskMap jumpProjector_;
  TaskMap fullTimeSlice_;
  TaskMap correctionSync_;
  TaskMap seedAssembler_;
  TaskMap leftSeedSync_;

  SeedMap mainSeed_; 

  DISALLOW_COPY_AND_ASSIGN(LocalNetwork);
};

} /* end namespace Hts */ } /* end namespace Pita */

#endif /* PITA_HTS_LOCALNETWORK_H */
