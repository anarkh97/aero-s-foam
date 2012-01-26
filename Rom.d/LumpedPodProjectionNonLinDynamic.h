#ifndef ROM_LUMPEDPODPROJECTIONNONLINDYNAMIC_H
#define ROM_LUMPEDPODPROJECTIONNONLINDYNAMIC_H

#include "PodProjectionNonLinDynamic.h"

#include <map>

namespace Rom {

class LumpedPodProjectionNonLinDynamic : public PodProjectionNonLinDynamic {
public:
  explicit LumpedPodProjectionNonLinDynamic(Domain *);

  virtual void preProcess();

private:
  // Overriden
  virtual void getStiffAndForceFromDomain(GeomState &geomState, Vector &elementInternalForce,
                                          Corotator **allCorot, FullSquareMatrix *kelArray,
                                          Vector &residual, double lambda, double time, GeomState *refState);
  
  std::map<int, double> packedElementWeights_;
  void buildPackedElementWeights();
};

} /* end namespace Rom */

#endif /* ROM_LUMPEDPODPROJECTIONNONLINDYNAMIC_H */
