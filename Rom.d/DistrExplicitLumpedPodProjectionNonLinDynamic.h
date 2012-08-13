#ifndef ROM_DISTREXPLICITLUMPEDPODPROJECTIONNONLINDYNAMIC_H
#define ROM_DISTREXPLICITLUMPEDPODPROJECTIONNONLINDYNAMIC_H

#include "DistrExplicitPodProjectionNonLinDynamicBase.h"

#include <vector>
#include <map>

namespace Rom {

class DistrExplicitLumpedPodProjectionNonLinDynamic : public DistrExplicitPodProjectionNonLinDynamicBase {
public:
  explicit DistrExplicitLumpedPodProjectionNonLinDynamic(Domain *);

  // Overriding via hiding
  void preProcess(); // Additional pre-processing
  void getInternalForce(DistrVector &d, DistrVector &f, double t); // Alternate internal force computation

private:
  void buildPackedElementWeights();
  
  void subGetWeightedInternalForceOnly(int iSub, DistrVector &f, double t);
  void subBuildPackedElementWeights(int iSub);

  std::vector<std::map<int, double> > packedElementWeights_;
};

} // end namespace Rom

#endif /* ROM_DISTREXPLICITLUMPEDPODPROJECTIONNONLINDYNAMIC_H */