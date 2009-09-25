#include "StaticSliceStrategy.h"

namespace Pita { namespace Hts {

TaskRank
StaticSliceStrategy::task(const SliceId & id) const {
  switch (id.type()) {
    case UNDEFINED_SLICE:
      return TaskRank(-1);
    case FORWARD_HALF_SLICE:  // Fall through
    case BACKWARD_HALF_SLICE: // Fall through
    case HEAD_FULL_SLICE:     // Fall through
    case TAIL_FULL_SLICE:
      return TaskRank(id.rank().value());
    default:
      break; // Fall through
  }
  throw Fwk::InternalException();
}

StaticSliceStrategy::SliceIdIterator
StaticSliceStrategy::slice(TaskRank workLoadId) const {
  std::vector<SliceId> slice;
  slice.reserve(4);
  slice.push_back(SliceId(FORWARD_HALF_SLICE, HalfSliceRank(workLoadId)));
  slice.push_back(SliceId(BACKWARD_HALF_SLICE, HalfSliceRank(workLoadId)));
  slice.push_back(SliceId(HEAD_FULL_SLICE, HalfSliceRank(workLoadId)));
  slice.push_back(SliceId(TAIL_FULL_SLICE, HalfSliceRank(workLoadId)));
  return iteratorNew(slice);
}

} /* end namespace Hts */ } /* end namespace Pita */
