#include "ReducedFullTimeSliceImpl.h"

#include <Math.d/FullSquareMatrix.h>

namespace Pita {

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
  /*if (iteration() == i) {
    return;
  }

  if (!correction() || !jump() || !nextCorrection()) {
    throw Exception("in  ReducedFullTimeSliceImpl::iterationIs -- Missing seeds");
  }

  if (jump()->iteration() != i || correction()->iteration() != i) {
    throw RangeException("in  ReducedFullTimeSliceImpl::iterationIs -- Iteration mismatch");
  }*/

  // result <- normalMatrix^{-1} * (reprojectionMatrix * correction + jump)
  Vector result = jump()->state();

  /*log() << "jump " << jump()->name() << "=\n";
  for (int i = 0; i < result.size(); ++i) {
    log() << result[i] << " ";
  }
  log() << "\n";*/

  /*log() << "corr " << correction()->name() << " =\n";
  for (int i = 0; i < correction()->state().size(); ++i) {
    log() << correction()->state()[i] << " ";
  }
  log() << "\n";*/
  
  //log() << "result.size() = " << result.size() << "\n";
  //log() << "correction()->state().size() = " << const_cast<Vector &>(correction()->state()).size() << "\n";
  //log() << "reprojectionMatrix()->size() = " << reprojectionMatrix()->dim() << "\n";

  //log() << "In RFTS " << toString(headHalfSlice()) << "\n";

  const_cast<FullSquareMatrix *>(reprojectionMatrix())->multiply(
      const_cast<Vector &>(correction()->state()),
      result,
      1.0,
      FullSquareMatrix::TRANSPOSED);
  
  /*log() << "rhs =\n";
  for (int i = 0; i < result.size(); ++i) {
    log() << result[i] << " ";
  }
  log() << "\n";*/

  projectionSolver()->solution(result);

  /*log() << "next_corr " << nextCorrection()->name() << "=\n";
  for (int i = 0; i < result.size(); ++i) {
    log() << result[i] << " ";
  }
  log() << "\n";*/

  ReducedSeed::Status cvgStatus =
    (correction()->status() == ReducedSeed::CONVERGED) && (jump()->status() == ReducedSeed::CONVERGED) ?
    ReducedSeed::CONVERGED :
    ReducedSeed::ACTIVE;

  nextCorrection()->stateIs(result);
  nextCorrection()->statusIs(cvgStatus);
  //nextCorrection()->iterationIs(ir);
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

} // end namespace Pita
