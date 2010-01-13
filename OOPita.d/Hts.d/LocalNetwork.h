#ifndef PITA_HTS_LOCALNETWORK_H
#define PITA_HTS_LOCALNETWORK_H

#include "Fwk.h"
#include "Types.h"

#include "SliceMapping.h"

#include "../Seed.h"
#include "../RemoteState.h"

#include "HalfTimeSlice.h"
#include "../JumpProjector.h"
#include "ReducedFullTimeSlice.h"
#include "../UpdatedSeedAssembler.h"
#include "CorrectionTimeSlice.h"

#include "../SeedErrorEvaluator.h"

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
  FullSliceCount totalFullSlices() const { return mapping_->totalFullSlices(); }
  CpuCount availableCpus() const { return mapping_->availableCpus(); }
  HalfSliceCount maxWorkload() const { return mapping_->maxWorkload(); }
  CpuRank localCpu() const { return localCpu_; }

  /* Current state */
  HalfSliceRank firstActiveSlice() const { return mapping_->firstActiveSlice(); }
  HalfSliceRank firstInactiveSlice() const { return mapping_->firstInactiveSlice(); }
  HalfSliceCount convergedSlices() const { return mapping_->convergedSlices(); }
  void convergedSlicesInc();

  /* Local network elements */
  TaskList halfTimeSlices() const;
  TaskList activeHalfTimeSlices() const;
  TaskList activeJumpProjectors() const;
  TaskList activeLeftSeedSyncs() const;
  TaskList activeFullTimeSlices() const;
  TaskList activeCorrectionSyncs() const;
  TaskList activeSeedAssemblers() const;

  TaskList activeCoarseTimeSlices() const;
  TaskList activeFullCorrectionSyncs() const;
  
  SeedMap mainSeeds() const;
  MainSeedMap activeMainSeeds() const;

  LocalNetwork(SliceMapping * mapping,
               CpuRank localCpu,
               HalfTimeSlice::Manager * sliceMgr,
               JumpProjector::Manager * jumpProjMgr,
               ReducedFullTimeSlice::Manager * ftsMgr,
               UpdatedSeedAssembler::Manager * usaMgr,
               CorrectionTimeSlice::Manager * ctsMgr, // Can be NULL
               RemoteState::Manager * commMgr,
               SeedErrorEvaluator::Manager * jumpErrorMgr);

protected:
  void init();

  Seed * getSeed(const SeedId & id);
  ReducedSeed * getReducedSeed(const SeedId & id);

private:
  SliceMapping::Ptr mapping_;
  CpuRank localCpu_;

  HalfTimeSlice::Manager::Ptr sliceMgr_;
  JumpProjector::Manager::Ptr jumpProjMgr_;
  ReducedFullTimeSlice::Manager::Ptr ftsMgr_;
  UpdatedSeedAssembler::Manager::Ptr usaMgr_;
  CorrectionTimeSlice::Manager::Ptr ctsMgr_;

  Seed::Manager::Ptr seedMgr_;
  ReducedSeed::Manager::Ptr reducedSeedMgr_;

  RemoteState::Manager::Ptr commMgr_;
  
  HTSList halfTimeSlice_;
  TaskMap jumpProjector_;
  TaskMap evenFullTimeSlice_;
  TaskMap oddFullTimeSlice_;
  TaskMap correctionSync_;
  TaskMap seedAssembler_;
  TaskMap leftSeedSync_;

  TaskMap evenCoarseTimeSlice_;
  TaskMap oddCoarseTimeSlice_;
  TaskMap fullCorrectionSync_;
  
  SeedMap mainSeed_; 

  SeedErrorEvaluator::Manager::Ptr jumpErrorMgr_;

  DISALLOW_COPY_AND_ASSIGN(LocalNetwork);
};

} /* end namespace Hts */ } /* end namespace Pita */

#endif /* PITA_HTS_LOCALNETWORK_H */
