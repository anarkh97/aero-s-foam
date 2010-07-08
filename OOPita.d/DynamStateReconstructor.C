#include "DynamStateReconstructor.h"

namespace Pita {

DynamStateReconstructor::DynamStateReconstructor(const DynamStateBasis * rb) :
  finalState_(),
  reconstructionBasis_(rb)
{}

void
DynamStateReconstructor::reducedBasisComponentsIs(const Vector & c) {
  int rbs = static_cast<int>(reducedBasisSize());
  if (c.size() != rbs) {
    throw Fwk::RangeException("Dimension mismatch");
  }
  DynamState result(vectorSize(), 0.0);
  for (int i = 0; i < rbs; ++i) {
    if (c[i] != 0.0) {
      result.linAdd(c[i], reconstructionBasis_->state(i));
    }
  }
  setFinalState(result);
}

void
DynamStateReconstructor::reconstructionBasisIs(const DynamStateBasis * rb) {
  if (rb != reconstructionBasis_.ptr()) {
    setFinalState(DynamState());
    reconstructionBasis_ = rb;
  }
}

} // end namespace Pita
