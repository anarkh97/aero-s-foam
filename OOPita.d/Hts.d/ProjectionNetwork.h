#ifndef PITA_HTS_PROJECTIONNETWORK_H
#define PITA_HTS_PROJECTIONNETWORK_H

#include "Fwk.h"
#include "Types.h"

#include "AffineBasisCollector.h"
#include "../DynamStateBasis.h"
#include "../RankDeficientSolver.h"
#include <Math.d/FullSquareMatrix.h>

namespace Pita { namespace Hts {

class ProjectionNetwork : public Fwk::PtrInterface<ProjectionNetwork> {
public:
  virtual size_t reducedBasisSize() const = 0;

  virtual AffineBasisCollector * collector() const = 0;

  virtual const DynamStateBasis * projectionBasis() const = 0;
  virtual const DynamStateBasis * propagatedBasis() const = 0;
  virtual const RankDeficientSolver * normalMatrixSolver() const = 0;
  virtual const FullSquareMatrix * reprojectionMatrix() const = 0;

  EXPORT_PTRINTERFACE_TYPES(ProjectionNetwork);

protected:
  ProjectionNetwork() {}

private:
  DISALLOW_COPY_AND_ASSIGN(ProjectionNetwork);
};

} /* end namespace Hts */ } /* end namespace Pita */

#endif /* PITA_HTS_PROJECTIONNETWORK_H */
