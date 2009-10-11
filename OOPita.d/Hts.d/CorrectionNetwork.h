#ifndef PITA_HTS_CORRECTIONNETWORK_H
#define PITA_HTS_CORRECTIONNETWORK_H

#include "Fwk.h"
#include "Types.h"

#include "BasisCollector.h"

#include "../JumpProjector.h"
#include "../UpdatedSeedAssembler.h"
#include "ReducedFullTimeSlice.h"

namespace Pita { namespace Hts {

class CorrectionNetwork : public Fwk::PtrInterface<CorrectionNetwork> {
public:
  virtual size_t reducedBasisSize() const = 0;

  virtual BasisCollector * collector() const = 0;

  virtual JumpProjector::Manager * jumpProjectorMgr() const = 0;
  virtual UpdatedSeedAssembler::Manager * updatedSeedAssemblerMgr() const = 0;
  virtual ReducedFullTimeSlice::Manager * fullTimeSliceMgr() const = 0;

  EXPORT_PTRINTERFACE_TYPES(CorrectionNetwork);

protected:
  CorrectionNetwork() {}

private:
  DISALLOW_COPY_AND_ASSIGN(CorrectionNetwork);
};

} /* end namespace Hts */ } /* end namespace Pita */

#endif /* PITA_HTS_CORRECTIONNETWORK_H */
