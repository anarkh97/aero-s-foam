#ifndef PITA_TIMESLICEIMPL_H
#define PITA_TIMESLICEIMPL_H

#include "TimeSlice.h"

#include <iostream>

namespace Pita {

class TimeSliceImpl : public TimeSlice {
public:
  typedef Fwk::Ptr<TimeSliceImpl> Ptr;
  typedef Fwk::Ptr<const TimeSliceImpl> PtrConst;

  // Built-in notification from InitialInterface
  virtual void onSeed() {}
  virtual void onStatus() {}
  
  class InitialInterface;
  class FinalInterface;

protected:
  TimeSliceImpl(SliceRank rank, Status status); 
};

} // end namespace Pita

#endif /* PITA_TIMESLICEIMPL_H */
