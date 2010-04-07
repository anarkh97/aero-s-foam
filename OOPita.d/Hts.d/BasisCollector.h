#ifndef PITA_HTS_BASISCOLLECTOR_H
#define PITA_HTS_BASISCOLLECTOR_H

#include "Fwk.h"
#include "Types.h"

#include "HalfTimeSlice.h"
#include "../IntegratorPropagator.h"

namespace Pita { namespace Hts {

class BasisCollector : public Fwk::PtrInterface<BasisCollector> {
public:
  EXPORT_PTRINTERFACE_TYPES(BasisCollector);
  
  virtual IntegratorPropagator * source(const HalfSliceId &) const = 0;
  virtual size_t sourceCount() const = 0;
  virtual IntegratorPropagator * sourceNew(const HalfSliceId &, DynamTimeIntegrator *) = 0;
  virtual void sourceDel(const HalfSliceId &) = 0;

protected:
  BasisCollector() {}

private:
  DISALLOW_COPY_AND_ASSIGN(BasisCollector);
};

} /* end namespace Hts */ } /* end namespace Pita */

#endif /* PITA_HTS_BASISCOLLECTOR_H */
