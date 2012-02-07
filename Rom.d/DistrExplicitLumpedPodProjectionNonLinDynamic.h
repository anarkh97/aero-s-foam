#ifndef ROM_DISTREXPLICITLUMPEDPODPROJECTIONNONLINDYNAMIC_H
#define ROM_DISTREXPLICITLUMPEDPODPROJECTIONNONLINDYNAMIC_H

#include <Paral.d/MDDynam.h>

#include "DistrVecBasis.h"

#include <vector>
#include <map>

namespace Rom {

class DistrExplicitLumpedPodProjectionNonLinDynamic : public MultiDomainDynam {
public:
  explicit DistrExplicitLumpedPodProjectionNonLinDynamic(Domain *);

  // Overriding via hiding
  void preProcess(); // Required additional pre-processing
  MDDynamMat * buildOps(double, double, double); // Modified matrix solver
  void getInternalForce(DistrVector &d, DistrVector &f, double t); // Alternate internal force computation

  ~DistrExplicitLumpedPodProjectionNonLinDynamic();

private:
  void buildPackedElementWeights();
  
  void subGetWeightedInternalForceOnly(int iSub, DistrVector &f, double t);
  void subBuildPackedElementWeights(int iSub);

  DistrVecBasis projectionBasis_;
  std::vector<std::map<int, double> > packedElementWeights_;

  // Disallow copy and assignment
  DistrExplicitLumpedPodProjectionNonLinDynamic(const DistrExplicitLumpedPodProjectionNonLinDynamic &);
  DistrExplicitLumpedPodProjectionNonLinDynamic &operator=(const DistrExplicitLumpedPodProjectionNonLinDynamic &);
};

} // end namespace Rom

#endif /* ROM_DISTREXPLICITLUMPEDPODPROJECTIONNONLINDYNAMIC_H */
