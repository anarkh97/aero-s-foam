#ifndef PITA_HTS_STATICSLICESTRATEGY_H
#define PITA_HTS_STATICSLICESTRATEGY_H

#include "SliceStrategy.h"

namespace Pita { namespace Hts {

class StaticSliceStrategy : public SliceStrategy {
public:
  EXPORT_PTRINTERFACE_TYPES(StaticSliceStrategy);

  virtual TaskRank task(const SliceId & id) const; // Overriden
  virtual SliceIdIterator slice(TaskRank workLoadId) const; // Overriden
  
  static Ptr New() {
    return new StaticSliceStrategy();
  }

protected:
  StaticSliceStrategy() {}
};

} /* end namespace Hts */ } /* end namespace Pita */

#endif /* PITA_HTS_STATICSLICESTRATEGY_H */
