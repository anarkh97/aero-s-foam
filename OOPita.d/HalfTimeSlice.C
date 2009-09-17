#include "HalfTimeSlice.h"

namespace Pita {

HalfTimeSlice::HalfTimeSlice(HalfSliceRank r, HalfTimeSlice::Direction d) :
  rank_(r),
  direction_(d),
  //phase_(),
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

} // end namespace Pita

