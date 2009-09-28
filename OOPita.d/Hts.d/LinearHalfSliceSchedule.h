#ifndef PITA_HTS_LINEARHALFSLICESCHEDULE_H
#define PITA_HTS_LINEARHALFSLICESCHEDULE_H

#include "HalfSliceSchedule.h"

namespace Pita { namespace Hts {

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

} /* end namespace Hts */ } /* end namespace Pita */

#endif /* PITA_HTS_LINEARHALFSLICESCHEDULE_H */
