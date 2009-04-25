#ifndef PITA_HALFSLICEBASISCOLLECTOR_H
#define PITA_HALFSLICEBASISCOLLECTOR_H

#include "Fwk.h"
#include "Types.h"

#include "HalfTimeSlice.h"
#include "IntegratorPropagator.h"

namespace Pita {

class HalfSliceBasisCollector : public Fwk::PtrInterface<HalfSliceBasisCollector> {
public:
  EXPORT_PTRINTERFACE_TYPES(HalfSliceBasisCollector);
  
  virtual IntegratorPropagator * source(const HalfSliceId &) const = 0;
  virtual size_t sourceCount() const = 0;
  virtual IntegratorPropagator * sourceNew(const HalfSliceId &, DynamTimeIntegrator *) = 0;
  virtual void sourceDel(const HalfSliceId &) = 0;

protected:
  HalfSliceBasisCollector() {}

private:
  DISALLOW_COPY_AND_ASSIGN(HalfSliceBasisCollector);
};

} // end namespace Pita

#endif /* PITA_HALFSLICEBASISCOLLECTOR_H */
