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
  void updateDisplacement(DistrVector& dinc, DistrVector& d_n);
  void getInternalForce(DistrVector &d, DistrVector &f, double t, int tIndex); // Alternate internal force computation
  MDDynamMat * buildOps(double, double, double);

private:
  void buildPackedElementWeights();
  void subUpdateWeightedNodesOnly(int iSub, DistrVector &v);
  void subGetWeightedInternalForceOnly(int iSub, DistrVector &f, double &t, int &tIndex);
  void subBuildPackedElementWeights(int iSub);

  std::vector<std::map<int, double> > packedElementWeights_;
  std::vector<std::vector<int> > packedWeightedNodes_;
};

} // end namespace Rom

#endif /* ROM_DISTREXPLICITLUMPEDPODPROJECTIONNONLINDYNAMIC_H */
