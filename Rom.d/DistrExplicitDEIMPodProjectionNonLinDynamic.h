#ifndef ROM_DISTREXPLICITDEIMPODPROJECTIONNONLINDYNAMIC_H
#define ROM_DISTREXPLICITDEIMPODPROJECTIONNONLINDYNAMIC_H

#include "DistrExplicitLumpedPodProjectionNonLinDynamic.h"

#include <vector>
#include <map>
#include <utility>

namespace Rom {

class DistrExplicitDEIMPodProjectionNonLinDynamic : public DistrExplicitLumpedPodProjectionNonLinDynamic {
public:
  explicit DistrExplicitDEIMPodProjectionNonLinDynamic(Domain *);

  // Overriding via hiding
  void preProcess(); // Additional pre-processing
  void getInternalForce(DistrVector &d, DistrVector &f, double t, int tIndex); // Alternate internal force computation

private:
  void buildInterpolationBasis();
  void subBuildInterpolationBasis(int iSub, std::vector< std::vector<std::pair<int,int> > > &maskedIndicesBuf);
  void subGetKtimesU(int isub, DistrVector &d, DistrVector &f);

  std::vector<std::map<int, double> > packedElementWeights_;
  std::vector<std::vector<int> > packedWeightedNodes_;

  DistrVecBasis deimBasis_;
  DistrVector * lin_fInt;
  GenFullSquareMatrix<double> **kelArrayCopy;
};

} // end namespace Rom

#endif /* ROM_DISTREXPLICITLUMPEDPODPROJECTIONNONLINDYNAMIC_H */
