#ifndef PITA_HTS_CORRECTIONNETWORK_H
#define PITA_HTS_CORRECTIONNETWORK_H

#include "Fwk.h"
#include "Types.h"

#include "BasisCollector.h"
#include "../DynamStateBasis.h"
#include "../RankDeficientSolver.h"
#include <Math.d/FullSquareMatrix.h>

namespace Pita { namespace Hts {

class CorrectionNetwork : public Fwk::PtrInterface<CorrectionNetwork> {
public:
  virtual size_t reducedBasisSize() const = 0;

  virtual BasisCollector * collector() const = 0;

  virtual const DynamStateBasis * projectionBasis() const = 0;
  virtual const DynamStateBasis * propagatedBasis() const = 0;
  virtual const RankDeficientSolver * normalMatrixSolver() const = 0;
  virtual const FullSquareMatrix * reprojectionMatrix() const = 0;

  EXPORT_PTRINTERFACE_TYPES(CorrectionNetwork);

protected:
  CorrectionNetwork() {}

private:
  DISALLOW_COPY_AND_ASSIGN(CorrectionNetwork);
};

} /* end namespace Hts */ } /* end namespace Pita */

#endif /* PITA_HTS_CORRECTIONNETWORK_H */
