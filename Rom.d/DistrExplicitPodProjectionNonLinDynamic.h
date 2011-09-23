#ifndef ROM_DISTREXPLICITPODPROJECTIONNONLINDYNAMIC_H
#define ROM_DISTREXPLICITPODPROJECTIONNONLINDYNAMIC_H

#include <Paral.d/MDDynam.h>

#include "DistrVecBasis.h"

namespace Rom {

template <typename Scalar> class GenDistrGalerkinProjectionSolver;

class DistrExplicitPodProjectionNonLinDynamic : public MultiDomainDynam {
public:
  explicit DistrExplicitPodProjectionNonLinDynamic(Domain *);

  // Overriding via hiding
  void preProcess(); // Required additional pre-processing
  MDDynamMat * buildOps(double, double, double); // Modified matrix solver

private:
  DistrVecBasis projectionBasis_;

  // Disallow copy and assignment
  DistrExplicitPodProjectionNonLinDynamic(const DistrExplicitPodProjectionNonLinDynamic &);
  DistrExplicitPodProjectionNonLinDynamic &operator=(const DistrExplicitPodProjectionNonLinDynamic &);
};

} // end namespace Rom

#endif /* ROM_DISTREXPLICITPODPROJECTIONNONLINDYNAMIC_H */
