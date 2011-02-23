#ifndef ROM_MESHSAMPLINGDRIVER_H
#define ROM_MESHSAMPLINGDRIVER_H

#include "DriverInterface.h"

#include "FileNameInfo.h"
#include "VecBasis.h"

#include <Driver.d/GeoSource.h> 
#include <Utils.d/DistHelper.h>

#include <vector>
#include <ostream>

class Domain;
class VecNodeDof6Conversion;
class NodalRestrictionMapping;
class SampledMeshRenumbering;
class MeshDesc;

RomDriverInterface *meshSamplingDriverNew(Domain *);

class MeshSamplingDriver : public RomDriverInterface {
public:
  virtual void solve();
  
  explicit MeshSamplingDriver(Domain *);

private:
  void preProcess();
  void buildDomainCdsa(); 

  // Input Pod 
  static void readPodOperators(const FileNameInfo &, const VecNodeDof6Conversion &, int maxDim, std::vector<VecBasis> &);
  static void readPodBasis(const BasisFileId &, const VecNodeDof6Conversion &, int maxDim, VecBasis &);

  // Reduced operators
  static void buildRestrictedProjection(const NodalRestrictionMapping &, const VecBasis &, VecBasis &);
  static void buildCombinedRestrictedProjection(const NodalRestrictionMapping &, const VecBasis &originPod, const VecBasis &targetPod, VecBasis &);

  // Output mesh
  static void outputMeshFile(const FileNameInfo &, const MeshDesc &);
  static std::string getMeshFilename(const FileNameInfo &);

  // Output Pod
  static void outputReducedOperators(const FileNameInfo &, const VecNodeDof6Conversion &, const NodalRestrictionMapping &, const VecBasis &jac, const VecBasis &res);
  static void outputExtendedPod(const BasisFileId &, const VecNodeDof6Conversion &, const NodalRestrictionMapping &, const VecBasis &);
  static void outputRestrictedPod(const BasisFileId &, const VecNodeDof6Conversion &, const NodalRestrictionMapping &, const VecBasis &);

  // Console output
  static void printReadPodInfo(BasisId::Type, int vectorCount);
  static void printGappyInfo(const std::vector<VecBasis> &podBases, const NodalRestrictionMapping &);

  Domain *domain_;
};

#endif /* ROM_MESHSAMPLINGDRIVER_H */
