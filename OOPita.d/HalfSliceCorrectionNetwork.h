#ifndef PITA_HALFSLICECORRECTIONNETWORK_H
#define PITA_HALFSLICECORRECTIONNETWORK_H

#include "Fwk.h"
#include "Types.h"

#include "HalfSliceBasisCollector.h"

#include "JumpProjector.h"
#include "UpdatedSeedAssembler.h"
#include "ReducedFullTimeSlice.h"

#include "DynamStateReductor.h"
#include "DynamStateReconstructor.h"

namespace Pita {

class HalfSliceCorrectionNetwork : public Fwk::PtrInterface<HalfSliceCorrectionNetwork> {
public:
  virtual size_t reducedBasisSize() const = 0;

  virtual HalfSliceBasisCollector * collector() const = 0;

  virtual JumpProjector::Manager * jumpProjectorMgr() const = 0;
  virtual UpdatedSeedAssembler::Manager * updatedSeedAssemblerMgr() const = 0;
  virtual ReducedFullTimeSlice::Manager * fullTimeSliceMgr() const = 0;

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
