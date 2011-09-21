#include "MeshSamplingDriver.h"

#include "RenumberingUtils.h"
#include "MeshDesc.h"
#include "VecBasisUtils.h"
#include "GreedyUtils.h"
#include "QrPseudoInversion.h"

#include "VecNodeDof6Map.h"
#include "VecNodeDof6Conversion.h"
#include "NodalRestrictionMapping.h"

#include "FileNameInfo.h"
#include "BasisFileStream.h"
#include "VecBasisFile.h"

#include "SimpleBuffer.h"

#include <Driver.d/Domain.h>

#include <vector>
#include <iterator>
#include <cstddef>
#include <fstream>
#include <cctype>
#include <algorithm>
#include <string>
#include <cassert>
#include <memory>

extern GeoSource *geoSource;

namespace Rom {

namespace { // anonymous

  template <typename T>
  inline
  const T *
  toArray(const std::vector<T> &v) {
    return !v.empty() ? &v[0] : NULL;
  }

} // end anonymous namespace

MeshSamplingDriver::MeshSamplingDriver(Domain *domain) :
  domain_(domain)
{}

void
MeshSamplingDriver::solve() {
  preProcess();
  
  // Read residual & jacobian POD bases 
  const FileNameInfo fileInfo;
  const VecNodeDof6Conversion vecDofConversion(*domain_->getCDSA());
  const int podSizeMax = domain_->solInfo().maxSizePodRom;
  std::vector<VecBasis> projectionPodBases;
  readPodOperators(fileInfo, vecDofConversion, podSizeMax, projectionPodBases);
  const VecBasis &residualPod = projectionPodBases[0];
  const VecBasis &jacobianPod = projectionPodBases[1];
  
  // Perform greedy selection of sample nodes
  VecNodeDof6Map vecDofMap(*domain_->getCDSA());
  std::vector<int> sampleNodeIds;
  greedy_sampling(projectionPodBases.begin(), projectionPodBases.end(),
                  geoSource->sampleNodeBegin(), geoSource->sampleNodeEnd(),
                  std::back_inserter(sampleNodeIds),
                  vecDofMap,
                  domain_->solInfo().aspectRatioPodRom);

  // Determine mapping between full and reduced mesh
  std::auto_ptr<Connectivity> elemToNode(new Connectivity(geoSource->getElemSet()));
  std::auto_ptr<Connectivity> nodeToElem(elemToNode->reverse());
  std::auto_ptr<Connectivity> nodeToNode(nodeToElem->transcon(elemToNode.get()));

  SampledMeshRenumbering meshRenumbering(sampleNodeIds.begin(), sampleNodeIds.end(),
                                         *nodeToNode, *nodeToElem);

  const NodalRestrictionMapping sampledMeshMapping(*domain_->getCDSA(),
                                                   meshRenumbering.sampleNodeIds().begin(),
                                                   meshRenumbering.sampleNodeIds().end());
  printGappyInfo(projectionPodBases, sampledMeshMapping);
 
  // Jacobian -> A
  VecBasis jacobianProjection;
  buildRestrictedProjection(sampledMeshMapping, jacobianPod, jacobianProjection);

  // (Residual, Jacobian) -> B
  VecBasis residualProjection;
  buildCombinedRestrictedProjection(sampledMeshMapping, residualPod, jacobianPod, residualProjection);

  // Assemble and output reduced mesh
  const MeshDesc reducedMesh(domain_, geoSource, meshRenumbering);
  outputMeshFile(fileInfo, reducedMesh);
 
  // Output sampled jacobian & residual POD bases
  const DofSetArray reducedDsa(meshRenumbering.reducedNodeIds().size(), const_cast<Elemset &>(reducedMesh.elements()));
  const ConstrainedDSA reducedCdsa(const_cast<DofSetArray &>(reducedDsa),
                                   reducedMesh.dirichletBConds().size(),
                                   const_cast<BCond *>(toArray(reducedMesh.dirichletBConds())));
  const VecNodeDof6Conversion reducedVecDofConversion(reducedCdsa);
  const NodalRestrictionMapping reducedSampledMeshMapping(reducedCdsa,
                                                          meshRenumbering.reducedSampleNodeIds().begin(),
                                                          meshRenumbering.reducedSampleNodeIds().end());
  outputReducedOperators(fileInfo, reducedVecDofConversion, reducedSampledMeshMapping, jacobianProjection, residualProjection);

  // Restrict state POD basis
  VecBasis statePod;
  readPodBasis(BasisFileId(fileInfo, BasisId::STATE, BasisId::POD), vecDofConversion, podSizeMax, statePod);
  const NodalRestrictionMapping reducedMeshMapping(*domain_->getCDSA(),
                                                   meshRenumbering.reducedNodeIds().begin(),
                                                   meshRenumbering.reducedNodeIds().end());
  outputRestrictedPod(BasisFileId(fileInfo, BasisId::STATE, BasisId::GAPPY_POD), reducedVecDofConversion, reducedMeshMapping, statePod);
}

void
MeshSamplingDriver::preProcess() {
  domain_->preProcessing();
  buildDomainCdsa();
}

void
MeshSamplingDriver::buildDomainCdsa() {
  const int numdof = domain_->numdof();
  SimpleBuffer<int> bc(numdof);
  SimpleBuffer<double> bcx(numdof);

  domain_->make_bc(bc.array(), bcx.array());
  domain_->make_constrainedDSA(bc.array());
}

void
MeshSamplingDriver::readPodOperators(const FileNameInfo &fileInfo,
                                     const VecNodeDof6Conversion &conversion,
                                     int maxDim,
                                     std::vector<VecBasis> &result) {
  const BasisId::Type types[] = { BasisId::RESIDUAL, BasisId::JACOBIAN };
  const int basisCount = sizeof(types) / sizeof(types[0]);

  result.resize(basisCount);
  for (std::vector<VecBasis>::iterator basisIt = result.begin();
      basisIt != result.end();
      ++basisIt) {
    const BasisId::Type t = types[std::distance(result.begin(), basisIt)];
    readPodBasis(BasisFileId(fileInfo, t, BasisId::POD), conversion, maxDim, *basisIt);
  }
}

void
MeshSamplingDriver::readPodBasis(const BasisFileId &basisId,
                                 const VecNodeDof6Conversion &conversion,
                                 int maxDim,
                                 VecBasis &result) {
  BasisInputStream in(basisId, conversion);
  const int vectorCount = maxDim ? std::min(maxDim, in.size()) : in.size();
  printReadPodInfo(basisId.type(), vectorCount);

  result.dimensionIs(vectorCount, in.vectorSize());
  readVectors(in, result.begin(), result.end());
}

void
MeshSamplingDriver::outputMeshFile(const FileNameInfo &fileInfo, const MeshDesc &mesh) {
  std::ofstream meshOut(getMeshFilename(fileInfo).c_str());
  meshOut << mesh;
}

std::string
MeshSamplingDriver::getMeshFilename(const FileNameInfo &fileInfo) {
  return fileInfo.prefix() + ".sampledmesh.inc";
}

void
MeshSamplingDriver::buildRestrictedProjection(const NodalRestrictionMapping &restriction,
                                              const VecBasis &originPod, VecBasis &result) {
  restrict_basis(restriction, originPod, result);
  QrPseudoInversion()(result);
}

void
MeshSamplingDriver::buildCombinedRestrictedProjection(const NodalRestrictionMapping &restriction,
                                                      const VecBasis &originPod, const VecBasis &targetPod,
                                                      VecBasis &result) {
  VecBasis preProjection;
  buildRestrictedProjection(restriction, originPod, preProjection);
  result.dimensionIs(targetPod.vectorCount(), preProjection.vectorSize());
  combine_projections(targetPod, originPod, preProjection, result);
}
  
void
MeshSamplingDriver::outputReducedOperators(const FileNameInfo &fileInfo,
                                           const VecNodeDof6Conversion &conversion,
                                           const NodalRestrictionMapping &restriction,
                                           const VecBasis &jacobianProjection,
                                           const VecBasis &residualProjection) {
  outputExtendedPod(BasisFileId(fileInfo, BasisId::JACOBIAN, BasisId::GAPPY_POD), conversion, restriction, jacobianProjection);
  outputExtendedPod(BasisFileId(fileInfo, BasisId::RESIDUAL, BasisId::GAPPY_POD), conversion, restriction, residualProjection);
}

void
MeshSamplingDriver::outputExtendedPod(const BasisFileId &fileId,
                                      const VecNodeDof6Conversion &conversion,
                                      const NodalRestrictionMapping &restriction,
                                      const VecBasis &basis) {
  BasisOutputStream out(fileId, conversion);
  writeExtendedVectors(out, basis, restriction);
}

void
MeshSamplingDriver::outputRestrictedPod(const BasisFileId &fileId,
                                        const VecNodeDof6Conversion &conversion,
                                        const NodalRestrictionMapping &restriction,
                                        const VecBasis &basis) {
  BasisOutputStream out(fileId, conversion);
  writeRestrictedVectors(out, basis, restriction);
}

void
MeshSamplingDriver::printReadPodInfo(BasisId::Type type, int vectorCount) {
  std::string name = toString(type);
  assert(!name.empty());
  name[0] = static_cast<char>(std::toupper(name[0]));
  filePrint(stderr, "%s POD subspace of dimension %d\n", name.c_str(), vectorCount);
}

void
MeshSamplingDriver::printGappyInfo(const std::vector<VecBasis> &podBases, const NodalRestrictionMapping &sampling) {
  filePrint(stderr, "Gappy mesh has %d sample nodes\n", sampling.sampleNodeCount());

  const int sampleLocationCount = sampling.restrictedInfo();
  const int podDimensionMax = vector_count_max(podBases.begin(), podBases.end());
  const double aspectRatio = static_cast<double>(sampleLocationCount) / podDimensionMax;
  filePrint(stderr, "Gappy aspect ratio = Sampled dofs / POD dimension = %d / %d = %g\n",
                    sampleLocationCount, podDimensionMax, aspectRatio);
}

} /* end namespace Rom */

Rom::DriverInterface *meshSamplingDriverNew(Domain *domain) {
  return new Rom::MeshSamplingDriver(domain);
}
