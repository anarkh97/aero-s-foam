#ifndef PITA_LINEARHALFSLICESCHEDULE_H
#define PITA_LINEARHALFSLICESCHEDULE_H

#include "HalfSliceSchedule.h"

namespace Pita {

class LinearHalfSliceSchedule : public HalfSliceSchedule {
public:
  typedef Fwk::Ptr<HalfSliceSchedule> Ptr;
  typedef Fwk::Ptr<const HalfSliceSchedule> PtrConst;

  static Ptr New(HalfSliceCount totalHalfSliceCount) {
    return new LinearHalfSliceSchedule(totalHalfSliceCount);
  }

protected:
  explicit LinearHalfSliceSchedule(HalfSliceCount nts);
};

} // end namespace Pita

#endif /* PITA_LINEARHALFSLICESCHEDULE_H */
