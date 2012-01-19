#ifndef ROM_DISTREXPLICITGAPPYNONLINDYNAMIC_H
#define ROM_DISTREXPLICITGAPPYNONLINDYNAMIC_H

#include <Paral.d/MDDynam.h>

#include "DistrVecBasis.h"

#include <memory>

namespace Rom {

template <typename Scalar> class GenDistrExplicitGappyProjectionSolver;

class DistrExplicitGappyNonLinDynamic : public MultiDomainDynam {
public:
  explicit DistrExplicitGappyNonLinDynamic(Domain *);

  // Overriding via hiding
  void preProcess(); // Required additional pre-processing
  MDDynamMat * buildOps(double, double, double); // Modified matrix solver

private:
  DistrVecBasis projectionBasis_;
  DistrVecBasis reconstructionBasis_;

  // Disallow copy and assignment
  DistrExplicitGappyNonLinDynamic(const DistrExplicitGappyNonLinDynamic &);
  DistrExplicitGappyNonLinDynamic &operator=(const DistrExplicitGappyNonLinDynamic &);
};

} // end namespace Rom

#endif /* ROM_DISTREXPLICITGAPPYNONLINDYNAMIC_H */
