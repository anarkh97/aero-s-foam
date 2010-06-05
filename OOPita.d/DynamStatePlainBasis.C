#include "DynamStatePlainBasis.h"

namespace Pita {

void
DynamStatePlainBasis::lastStateIs(const DynamState & ds) {
  if (ds.vectorSize() != this->vectorSize())
    throw Fwk::RangeException(); 
  addState(ds);
}

void
DynamStatePlainBasis::lastStateBasisIs(const DynamStateBasis * dsb) {
  if (dsb->vectorSize() != this->vectorSize())
    throw Fwk::RangeException();
  size_t statesToAdd = dsb->stateCount();
  for (size_t i = 0; i < statesToAdd; ++i)
    addState(dsb->state(i));
}

void
DynamStatePlainBasis::lastStateBasisIs(const DynamStatePlainBasis * dsb) {
  if (dsb->vectorSize() != this->vectorSize())
    throw Fwk::RangeException();
  addStateBasis(dsb);
}


} // end namespace Pita
