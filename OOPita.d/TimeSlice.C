#include "TimeSlice.h"

namespace Pita {

TimeSlice::TimeSlice(SliceRank rank, Status status) : 
  rank_(rank),
  status_(status),
  initialInterface_(NULL), 
  finalInterface_(NULL),
  notifiee_(NULL)
{}

} // end namespace Pita
