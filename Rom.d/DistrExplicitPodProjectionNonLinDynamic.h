#ifndef ROM_DISTREXPLICITPODPROJECTIONNONLINDYNAMIC_H
#define ROM_DISTREXPLICITPODPROJECTIONNONLINDYNAMIC_H

#include <Paral.d/MDDynam.h>

#include "DistrVecBasis.h"

#include <Problems.d/DynamProbTraits.h>

namespace Rom {

template <typename Scalar> class GenDistrGalerkinProjectionSolver;

class DistrExplicitPodProjectionNonLinDynamic : public MultiDomainDynam {
public:
  explicit DistrExplicitPodProjectionNonLinDynamic(Domain *);

  // Overriding via hiding
  void preProcess(); // Required additional pre-processing
  MDDynamMat * buildOps(double, double, double); // Modified matrix solver

  // Added functionality
  void forceSnapshotAdd(const DistrVector &s);

  ~DistrExplicitPodProjectionNonLinDynamic();

private:
  DistrVecBasis projectionBasis_;

  class SnapshotHandler;
  friend class SnapshotHandler;
  SnapshotHandler * snapshotHandler_;

  // Disallow copy and assignment
  DistrExplicitPodProjectionNonLinDynamic(const DistrExplicitPodProjectionNonLinDynamic &);
  DistrExplicitPodProjectionNonLinDynamic &operator=(const DistrExplicitPodProjectionNonLinDynamic &);
};

} // end namespace Rom

inline
void
handleForce(Rom::DistrExplicitPodProjectionNonLinDynamic &probDesc, DistrVector &f) {
  probDesc.forceSnapshotAdd(f);
}


#endif /* ROM_DISTREXPLICITPODPROJECTIONNONLINDYNAMIC_H */
