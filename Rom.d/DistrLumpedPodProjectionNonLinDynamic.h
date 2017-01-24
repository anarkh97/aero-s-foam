#ifndef ROM_DISTRLUMPEDPODPROJECTIONNONLINDYNAMIC_H
#define ROM_DISTRLUMPEDPODPROJECTIONNONLINDYNAMIC_H

#include "DistrPodProjectionNonLinDynamic.h"

#include <map>

namespace Rom {

class DistrLumpedPodProjectionNonLinDynamic : public DistrPodProjectionNonLinDynamic {
public:
  explicit DistrLumpedPodProjectionNonLinDynamic(Domain *);

  virtual void preProcess();
  virtual void updateStates(DistrModalGeomState *refState, DistrModalGeomState& geomState, double time);
  virtual void setLocalReducedMesh(int j);

private:
  void subBuildPackedElementWeights(int iSub);
  void subZeroStiffMats(int iSub);  
  void subUpdateStates(int iSub, DistrGeomState *refState, DistrGeomState *geomState, double time);
  double getStiffAndForce(DistrModalGeomState& geomState, DistrVector& residual, DistrVector& elementInternalForce,
                            double midtime=-1, DistrModalGeomState *refState = NULL, bool forceOnly = false);
  void subGetStiffAndForce(int iSub, DistrGeomState &geomState, DistrVector &res,
                             DistrVector &elemIntForce, double t, DistrGeomState *refState, bool forceOnly);

protected:
  std::vector<std::vector<std::map<int, double> > > packedElementWeights_;
  std::vector<std::vector<std::vector<int> > > localPackedWeightedNodes_;
  std::vector<std::vector<int> > packedWeightedNodes_;
  std::set<int> packedWeightedElems_;
  void buildPackedElementWeights();
  int localReducedMeshId_;
};

} /* end namespace Rom */

#endif /* ROM_LUMPEDPODPROJECTIONNONLINDYNAMIC_H */
