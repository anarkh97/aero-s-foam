#include "ReducedFullTimeSlice.h"

namespace Pita {

// ReducedFullTimeSliceHead implemetation

ReducedFullTimeSliceHead::ReducedFullTimeSliceHead(HalfSliceRank headRank, CpuRank tailCpu) :
  headHalfSlice_(headRank),
  tailHalfSlice_(HalfSliceRank(headRank.value() + 1)),
  tailCpu_(tailCpu),
  phase_(),
  rightPropagatedSeed_(NULL),
  leftPropagatedSeed_(NULL),
  updatedSeed_(NULL),
  correction_(NULL)
{}

// ReducedFullTimeSliceTail implemetation

ReducedFullTimeSliceTail::ReducedFullTimeSliceTail(HalfSliceRank tailRank, CpuRank headCpu) :
  headHalfSlice_(HalfSliceRank(tailRank.value() - 1)),
  tailHalfSlice_(tailRank),
  headCpu_(headCpu),
  nextLeftPropagatedSeed_(NULL),
  nextUpdatedSeed_(NULL),
  nextCorrection_(NULL)
{}

} // end namespace Pita
