#ifndef PITA_HTS_BASISCOLLECTOR_H
#define PITA_HTS_BASISCOLLECTOR_H

#include "Fwk.h"
#include "Types.h"

#include "../DynamPropagator.h"

#include "HalfTimeSlice.h"

namespace Pita { namespace Hts {

class BasisCollector : public Fwk::PtrInterface<BasisCollector> {
public:
  EXPORT_PTRINTERFACE_TYPES(BasisCollector);
  
  virtual DynamPropagator * source(const HalfSliceId & id) const = 0;
  virtual size_t sourceCount() const = 0;
  virtual void sourceIs(const HalfSliceId & id, DynamPropagator * source) = 0;

protected:
  BasisCollector() {}

private:
  DISALLOW_COPY_AND_ASSIGN(BasisCollector);
};

} /* end namespace Hts */ } /* end namespace Pita */

#endif /* PITA_HTS_BASISCOLLECTOR_H */
