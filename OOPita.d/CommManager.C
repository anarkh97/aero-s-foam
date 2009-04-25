#include "CommManager.h"

namespace Pita {

CommManager::CommManager(TimeSliceMapping * sliceMap, CpuRank myCpuRank, CpuRank nextCpuRank, CpuRank previousCpuRank) :
  status_(available),
  sliceMap_(sliceMap),
  localCpuRank_(myCpuRank),
  nextCpuRank_(nextCpuRank),
  previousCpuRank_(previousCpuRank),
  notifiee_(NULL),
  localBroadcastBasis_(NULL),
  receivedBroadcastBasis_(NULL)
{}
 
} // end namespace Pita
