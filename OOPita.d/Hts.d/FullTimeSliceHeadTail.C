#include "FullTimeSliceHeadTail.h"

namespace Pita { namespace Hts {

// FullTimeSliceHead implemetation

FullTimeSliceHead::FullTimeSliceHead(HalfSliceRank headRank, CpuRank tailCpu) :
  headHalfSlice_(headRank),
  tailHalfSlice_(HalfSliceRank(headRank.value() + 1)),
  tailCpu_(tailCpu),
  phase_(),
  updatedSeed_(NULL),
  rightPropagatedSeed_(NULL)
{}

// FullTimeSliceTail implemetation

FullTimeSliceTail::FullTimeSliceTail(HalfSliceRank tailRank, CpuRank headCpu) :
  headHalfSlice_(HalfSliceRank(tailRank.value() - 1)),
  tailHalfSlice_(tailRank),
  headCpu_(headCpu),
  nextLeftPropagatedSeed_(NULL),
  nextUpdatedSeed_(NULL)
{}

} /* end namespace Hts */ } /* end namespace Pita */
