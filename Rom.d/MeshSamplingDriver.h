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
  static void readPodOperators(const FileNameInfo &,
                               const VecNodeDof6Conversion &,
                               std::vector<VecBasis> &, int podSizeMax);
  static void readPodBasis(VecBasis &, BasisId::Type, const FileNameInfo &,
                           const VecNodeDof6Conversion &, int maxDim);

  // Reduced operators
  static void buildRestrictedProjection(const NodalRestrictionMapping &, const VecBasis &, VecBasis &);
  static void buildCombinedRestrictedProjection(const NodalRestrictionMapping &,
                                                const VecBasis &originPod, const VecBasis &targetPod,
                                                VecBasis &);

  // Output mesh
  static void outputMeshFile(const FileNameInfo &, const MeshDesc &);
  static std::string getMeshFilename(const FileNameInfo &);

  // Output Pod
  static void outputReducedOperators(const NodalRestrictionMapping &, const FileNameInfo &, const VecNodeDof6Conversion &,
                                     const VecBasis &jac, const VecBasis &res);
  static void outputSampledPod(const VecBasis &, BasisId::Type, const FileNameInfo &,
                               const VecNodeDof6Conversion &, const NodalRestrictionMapping &);
  static void outputStatePodRestriction(const NodalRestrictionMapping &restriction,
                                        const FileNameInfo &fileInfo,
                                        const VecNodeDof6Conversion &conversion,
                                        const VecBasis &basis);

  // Console output
  static void printReadPodInfo(BasisId::Type, int vectorCount);
  static void printGappyInfo(const std::vector<VecBasis> &podBases, const NodalRestrictionMapping &);

  Domain *domain_;
};

#endif /* ROM_MESHSAMPLINGDRIVER_H */
