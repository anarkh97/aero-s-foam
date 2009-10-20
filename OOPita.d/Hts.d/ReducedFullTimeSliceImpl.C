#include "ReducedFullTimeSliceImpl.h"

#include <Math.d/FullSquareMatrix.h>

#include <cassert>

namespace Pita { namespace Hts {

// ReducedFullTimeSliceImpl implementation

ReducedFullTimeSliceImpl::ReducedFullTimeSliceImpl(HalfSliceRank headRank,
                                                   const FullSquareMatrix * reprojectionMatrix,
                                                   const RankDeficientSolver * projectionSolver):
  ReducedFullTimeSlice(headRank),
  reprojectionMatrix_(reprojectionMatrix),
  projectionSolver_(projectionSolver)
{}

void
ReducedFullTimeSliceImpl::iterationIs(IterationRank ir) {
  assert(jump()->iteration() == ir);
  assert(correction()->iteration() == ir);

  Vector result = jump()->state();
  
  // result == reprojectionMatrix * correction + jump 
  const_cast<FullSquareMatrix *>(reprojectionMatrix())->multiply(
      const_cast<Vector &>(correction()->state()),
      result,       
      1.0,
      FullSquareMatrix::TRANSPOSED);
  
  // result == normalMatrix^{-1} * (reprojectionMatrix * correction + jump)
  projectionSolver()->solution(result);

  ReducedSeed::Status cvgStatus =
    (correction()->status() == ReducedSeed::CONVERGED) && (jump()->status() == ReducedSeed::CONVERGED) ?
    ReducedSeed::CONVERGED :
    ReducedSeed::ACTIVE;

  nextCorrection()->stateIs(result);
  nextCorrection()->statusIs(cvgStatus);
  nextCorrection()->iterationIs(ir);
}

// ReducedFullTimeSliceImpl::Manager implementation

ReducedFullTimeSliceImpl::Manager::Manager(const FullSquareMatrix * defaultReprojectionMatrix,
                                           const RankDeficientSolver * defaultProjectionSolver) :
  defaultReprojectionMatrix_(defaultReprojectionMatrix),
  defaultProjectionSolver_(defaultProjectionSolver)
{}

ReducedFullTimeSliceImpl *
ReducedFullTimeSliceImpl::Manager::createNewInstance(const HalfSliceRank & r) {
  return new ReducedFullTimeSliceImpl(r, defaultReprojectionMatrix(), defaultProjectionSolver());
}

} /* end namespace Hts */ } /* end namespace Pita */
