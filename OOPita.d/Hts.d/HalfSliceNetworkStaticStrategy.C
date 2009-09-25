#include "OOPita.d/HalfSliceNetworkStaticStrategy.h"

namespace Pita {

HalfSliceRank
HalfSliceNetworkStaticStrategy::workLoadShare(const Hts::SliceId & id) const {
  switch (id.type()) {
    case Hts::UNDEFINED_SLICE:
      return HalfSliceRank(-1);
    case Hts::FORWARD_HALF_SLICE:  // Fall through
    case Hts::BACKWARD_HALF_SLICE: // Fall through
    case Hts::HEAD_FULL_SLICE:     // Fall through
    case Hts::TAIL_FULL_SLICE:
      return id.rank();
    default: 
      throw Fwk::InternalException();
  }
}

HalfSliceNetworkStrategy::SliceIteratorConst
HalfSliceNetworkStaticStrategy::slice(HalfSliceRank workLoadId) const {
  std::vector<Hts::SliceId> slice;
  slice.reserve(4);
  slice.push_back(Hts::SliceId(Hts::FORWARD_HALF_SLICE, workLoadId));
  slice.push_back(Hts::SliceId(Hts::BACKWARD_HALF_SLICE, workLoadId));
  slice.push_back(Hts::SliceId(Hts::HEAD_FULL_SLICE, workLoadId));
  slice.push_back(Hts::SliceId(Hts::TAIL_FULL_SLICE, workLoadId));
  return createIterator(slice);
}

} // end namespace Pita
