#ifndef PITA_HALFSLICECORRECTIONNETWORK_H
#define PITA_HALFSLICECORRECTIONNETWORK_H

#include "Fwk.h"
#include "Types.h"

#include "HalfSliceBasisCollector.h"
#include "DynamStateReductor.h"
#include "DynamStateReconstructor.h"

namespace Pita {

class HalfSliceCorrectionNetwork : public Fwk::PtrInterface<HalfSliceCorrectionNetwork> {
public:
  virtual HalfSliceBasisCollector * collector() const = 0;
  virtual DynamStateReductor::Manager * reductorMgr() const = 0;
  virtual DynamStateReconstructor::Manager * reconstructorMgr() const = 0;

  EXPORT_PTRINTERFACE_TYPES(HalfSliceCorrectionNetwork);

protected:
  HalfSliceCorrectionNetwork() {}

private:
  DISALLOW_COPY_AND_ASSIGN(HalfSliceCorrectionNetwork);
};

} // end namespace Pita

#endif /* PITA_HALFSLICECORRECTIONNETWORK_H */
