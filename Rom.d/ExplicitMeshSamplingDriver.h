#ifndef ROM_EXPLICITMESHSAMPLINGDRIVER_H
#define ROM_EXPLICITMESHSAMPLINGDRIVER_H

#include "DriverInterface.h"

#include "FileNameInfo.h"
#include "VecBasis.h"

#include <Driver.d/GeoSource.h> 
#include <Utils.d/DistHelper.h>

#include <vector>
#include <ostream>

class Domain;
template <typename Scalar> class GenSparseMatrix;

namespace Rom {

class VecNodeDof6Conversion;
class NodalRestrictionMapping;
class SampledMeshRenumbering;
class MeshDesc;

class ExplicitMeshSamplingDriver : public DriverInterface {
public:
  virtual void solve();
  
  explicit ExplicitMeshSamplingDriver(Domain *);

private:
  void preProcess();
  void buildDomainCdsa();
  GenSparseMatrix<double> * massMatrixNew();

  // Input Pod 
  static void readPodBasis(const BasisFileId &, const VecNodeDof6Conversion &, int maxDim, VecBasis &);

  // Reduced operators
  static void buildRestrictedProjection(const NodalRestrictionMapping &, const VecBasis &, VecBasis &);
  static void buildCombinedRestrictedProjection(const NodalRestrictionMapping &, const VecBasis &originPod, const VecBasis &targetPod, VecBasis &);

  // Output mesh
  static void outputMeshFile(const FileNameInfo &, const MeshDesc &);
  static std::string getMeshFilename(const FileNameInfo &);

  // Output Pod
  static void outputReducedOperators(const FileNameInfo &, const VecNodeDof6Conversion &, const NodalRestrictionMapping &, const VecBasis &force);
  static void outputExtendedPod(const BasisFileId &, const VecNodeDof6Conversion &, const NodalRestrictionMapping &, const VecBasis &);
  static void outputRestrictedPod(const BasisFileId &, const VecNodeDof6Conversion &, const NodalRestrictionMapping &, const VecBasis &);

  // Console output
  static void printReadPodInfo(BasisId::Type, int vectorCount);
  static void printGappyInfo(const VecBasis &podBasis, const NodalRestrictionMapping &);

  Domain *domain_;
};

} /* end namespace Rom */

Rom::DriverInterface *explicitMeshSamplingDriverNew(Domain *);

#endif /* ROM_EXPLICITMESHSAMPLINGDRIVER_H */
