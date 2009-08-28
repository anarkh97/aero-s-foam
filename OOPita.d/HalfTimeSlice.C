#include "HalfTimeSlice.h"

namespace Pita {

HalfTimeSlice::HalfTimeSlice(HalfSliceRank r, HalfTimeSlice::Direction d) :
  rank_(r),
  direction_(d),
  //phase_(),
  seed_(NULL),
  propagatedSeed_(NULL)
{}

OStream & operator<<(OStream & out, const HalfSliceId & id) {
  out << id.rank();
  out << (id.direction() == HalfTimeSlice::FORWARD ? 'F' : 'B');
  return out; 
}

} // end namespace Pita

