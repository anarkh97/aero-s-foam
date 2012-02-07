#ifndef ROM_DISTREXPLICITPODPROJECTIONNONLINDYNAMICBASE_H
#define ROM_DISTREXPLICITPODPROJECTIONNONLINDYNAMICBASE_H

#include <Paral.d/MDDynam.h>

#include "DistrVecBasis.h"

namespace Rom {

class DistrExplicitPodProjectionNonLinDynamicBase : public MultiDomainDynam {
public:
  explicit DistrExplicitPodProjectionNonLinDynamicBase(Domain *);

  // Overriding via hiding
  void preProcess(); // Additional pre-processing
  MDDynamMat * buildOps(double, double, double); // Reduced-order matrix solver

private:
  DistrVecBasis projectionBasis_;
  
  // Disallow copy and assignment
  DistrExplicitPodProjectionNonLinDynamicBase(const DistrExplicitPodProjectionNonLinDynamicBase &);
  DistrExplicitPodProjectionNonLinDynamicBase &operator=(const DistrExplicitPodProjectionNonLinDynamicBase &);
};

} // end namespace Rom

#endif /* ROM_DISTREXPLICITPODPROJECTIONNONLINDYNAMICBASE_H */
