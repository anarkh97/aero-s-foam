#ifndef ROM_UDEIM_SAMPLINGDRIVER_H
#define ROM_UDEIM_SAMPLINGDRIVER_H

#include "DriverInterface.h"

#include "Problems.d/DynamDescr.h"
#include "VecBasisOps.h"
#include "FileNameInfo.h"
#include "VecBasisFile.h"

#include "VecNodeDof6Conversion.h"

#ifdef USE_EIGEN3
#include <Eigen/LU>
#endif

namespace Rom {

class UDEIMSamplingDriver : public SingleDomainDynamic, public DriverInterface {
public:
  virtual void solve();
  
  explicit UDEIMSamplingDriver(Domain *);

private:

  void readInBasis(VecBasis &podBasis, BasisId::Type type, BasisId::Level level, int podSizeMax, bool normalized = false);
  template <typename Scalar> void writeBasisToFile(const VecBasis &OutputBasis, std::vector<Scalar> singularValue, BasisId::Type type,BasisId::Level);

  int  snapSize(BasisId::Type type, std::vector<int> &snapshotCounts);
  void readAndProjectSnapshots(BasisId::Type type, const int vectorSize, VecBasis &podBasis,
                          const VecNodeDof6Conversion *vecDofConversion,
                          std::vector<int> &snapshotCounts, std::vector<double> &timeStamps, VecBasis &config);
  void buildForceArray(VecBasis &forceBasis, const VecBasis &displac,
                       const VecBasis *veloc, const VecBasis *accel,std::vector<double> timeStamps_,
                       std::vector<int> snapshotCounts_);
  void OrthoForceSnap(VecBasis &forceBasis,std::vector<double> &SVs);  

  int  elementCount() const; 

  void writeProjForceSnap(); 

  void computeInterpIndices(VecBasis &forceBasis, std::vector<int> &maskIndices); 
  void computeAndWriteUDEIMBasis(VecBasis &forceBasis, std::vector<int> &maskIndices);  
  void writeSampledMesh(std::vector<int> &maskIndices);

  VecNodeDof6Conversion *converter;

  VecBasis podBasis_;

};

} /* end namespace Rom */

Rom::DriverInterface *udeimSamplingDriverNew(Domain *);

#endif /* ROM_BASIS_ORTHODRIVER_H */
