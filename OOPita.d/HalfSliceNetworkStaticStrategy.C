#include "HalfSliceNetworkStaticStrategy.h"

namespace Pita {

HalfSliceRank
HalfSliceNetworkStaticStrategy::workLoadShare(const Hs::SliceId & id) const {
  switch (id.type()) {
    case Hs::UNDEFINED_SLICE:
      return HalfSliceRank(-1);
    case Hs::FORWARD_HALF_SLICE:  // Fall through
    case Hs::BACKWARD_HALF_SLICE: // Fall through
    case Hs::HEAD_FULL_SLICE:     // Fall through
    case Hs::TAIL_FULL_SLICE:
      return id.rank();
    default: 
      throw Fwk::InternalException();
  }
}

HalfSliceNetworkStrategy::SliceIteratorConst
HalfSliceNetworkStaticStrategy::slice(HalfSliceRank workLoadId) const {
  std::vector<Hs::SliceId> slice;
  slice.reserve(4);
  slice.push_back(Hs::SliceId(Hs::FORWARD_HALF_SLICE, workLoadId));
  slice.push_back(Hs::SliceId(Hs::BACKWARD_HALF_SLICE, workLoadId));
  slice.push_back(Hs::SliceId(Hs::HEAD_FULL_SLICE, workLoadId));
  slice.push_back(Hs::SliceId(Hs::TAIL_FULL_SLICE, workLoadId));
  return createIterator(slice);
}

} // end namespace Pita
