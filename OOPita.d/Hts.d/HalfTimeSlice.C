#include "HalfTimeSlice.h"

namespace Pita { namespace Hts {

HalfTimeSlice::HalfTimeSlice(HalfSliceRank r, HalfTimeSlice::Direction d) :
  NamedTask(String("HalfTimeSlice ") + toString(HalfSliceId(r, d))),
  rank_(r),
  direction_(d),
  seed_(NULL),
  propagatedSeed_(NULL)
{}

OStream & operator<<(OStream & out, const HalfTimeSlice::Direction & d) {
  out << (d == HalfTimeSlice::FORWARD ? 'F' : 'B');
  return out;
}

OStream & operator<<(OStream & out, const HalfSliceId & id) {
  out << id.rank() << id.direction();
  return out; 
}

} /* end namespace Hts */ } /* end namespace Pita */

