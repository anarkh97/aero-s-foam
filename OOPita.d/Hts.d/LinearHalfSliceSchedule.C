#include "LinearHalfSliceSchedule.h"

namespace Pita { namespace Hts {

LinearHalfSliceSchedule::LinearHalfSliceSchedule(HalfSliceCount nts) :
  HalfSliceSchedule(PhaseRank(5), // localPropagation,
                    PhaseRank(3), // correctionPhase,
                    PhaseRank(4), // mainSeedSynchronization,
                    PhaseRank(1), // propagatedSeedSynchronization,
                    PhaseRank(2), // baseBuilding
                    nts)
{}

} /* end namespace Hts */ } /* end namespace Pita */
