#ifndef PITA_HALFSLICESCHEDULE_H
#define PITA_HALFSLICESCHEDULE_H

#include "Fwk.h"
#include "Types.h"
#include "../Activity.h"

namespace Pita { namespace Hts {

class HalfSliceSchedule : public Fwk::PtrInterface<HalfSliceSchedule> {
public:
  EXPORT_PTRINTERFACE_TYPES(HalfSliceSchedule);

  PhaseRank newIteration(HalfSliceRank hs = HalfSliceRank(-1)) const { return getAdjustedPhase(PhaseRank(0), hs); }
  PhaseRank localPropagation(HalfSliceRank hs = HalfSliceRank(-1)) const { return getAdjustedPhase(localPropagation_, hs); }
  PhaseRank correction(HalfSliceRank hs = HalfSliceRank(-1)) const { return getAdjustedPhase(correction_, hs); } 
  PhaseRank propagatedSeedSynchronization(HalfSliceRank hs = HalfSliceRank(-1)) const { return getAdjustedPhase(propagatedSeedSynchronization_, hs); }
  PhaseRank mainSeedSynchronization(HalfSliceRank hs = HalfSliceRank(-1)) const { return getAdjustedPhase(mainSeedSynchronization_, hs); }
  PhaseRank baseBuilding(HalfSliceRank hs = HalfSliceRank(-1)) const { return getAdjustedPhase(baseBuilding_, hs); }
  PhaseRank endIteration() const { return endIteration_; }

protected:
  HalfSliceSchedule(PhaseRank lp, PhaseRank c, PhaseRank mss, PhaseRank pss, PhaseRank bb, HalfSliceCount nts) :
    localPropagation_(lp),
    correction_(c),
    mainSeedSynchronization_(mss),
    propagatedSeedSynchronization_(pss),
    baseBuilding_(bb),
    endIteration_(),
    totalHalfSliceCount_(nts)
  {
    PhaseRank lastPhase = std::max(lp, std::max(c, std::max(mss, std::max(pss, bb))));
    endIteration_ = getAdjustedPhase(lastPhase, HalfSliceRank(0) + nts);
  }

  PhaseRank getAdjustedPhase(PhaseRank baseRank, HalfSliceRank hs) const {
    return PhaseRank(baseRank.value() * (totalHalfSliceCount_.value() + 1) + hs.value() + 1); 
  }

private:
  PhaseRank localPropagation_;
  PhaseRank correction_;
  PhaseRank mainSeedSynchronization_;
  PhaseRank propagatedSeedSynchronization_;
  PhaseRank baseBuilding_;
  PhaseRank endIteration_;
  HalfSliceCount totalHalfSliceCount_;

  DISALLOW_COPY_AND_ASSIGN(HalfSliceSchedule);
};

typedef HalfSliceSchedule Schedule;

} /* end namespace Hts */ } /* namespace Pita */

#endif /* PITA_HALFSLICESCHEDULE_H */
