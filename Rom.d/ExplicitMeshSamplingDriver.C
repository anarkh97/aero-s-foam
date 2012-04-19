#include "ExplicitMeshSamplingDriver.h"

#include "RenumberingUtils.h"
#include "MeshDesc.h"
#include "VecBasisUtils.h"
#include "VecBasisOps.h"
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
#include <memory>
#include <cassert>

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

ExplicitMeshSamplingDriver::ExplicitMeshSamplingDriver(Domain *domain) :
  domain_(domain)
{}

void
ExplicitMeshSamplingDriver::solve() {
  preProcess();
  
  // Read force POD basis
  const FileNameInfo fileInfo;
  const VecNodeDof6Conversion vecDofConversion(*domain_->getCDSA());
  const int podSizeMax = domain_->solInfo().maxSizePodRom;
  VecBasis forcePod[1];
  const BasisFileId forcePodId(fileInfo, BasisId::FORCE, BasisId::POD);
  readPodBasis(forcePodId, vecDofConversion, podSizeMax, forcePod[0]);
  
  // Perform greedy selection of sample nodes
  VecNodeDof6Map vecDofMap(*domain_->getCDSA());
  std::vector<int> sampleNodeIds;
  greedy_sampling(&forcePod[0], &forcePod[0] + 1,
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
  printGappyInfo(forcePod[0], sampledMeshMapping);
 
  // Normalize state POD basis
  VecBasis statePod;
  readPodBasis(BasisFileId(fileInfo, BasisId::STATE, BasisId::POD), vecDofConversion, podSizeMax, statePod);

  VecBasis normalizedStatePod;
  std::auto_ptr<const SparseMatrix> metric(massMatrixNew());
  renormalized_basis(*metric, statePod, normalizedStatePod);
  
  VecBasis forceProjection;
  buildCombinedRestrictedProjection(sampledMeshMapping, forcePod[0], normalizedStatePod, forceProjection);

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
  outputReducedOperators(fileInfo, reducedVecDofConversion, reducedSampledMeshMapping, forceProjection);

  const NodalRestrictionMapping reducedMeshMapping(*domain_->getCDSA(),
                                                   meshRenumbering.reducedNodeIds().begin(),
                                                   meshRenumbering.reducedNodeIds().end());
  outputRestrictedPod(BasisFileId(fileInfo, BasisId::STATE, BasisId::GAPPY_POD), reducedVecDofConversion, reducedMeshMapping, normalizedStatePod);
}

void
ExplicitMeshSamplingDriver::preProcess() {
  geoSource->setMRatio(0.0); // Enforce mass matrix lumping
  domain_->preProcessing();
  buildDomainCdsa();
  domain_->makeAllDOFs();
}

void
ExplicitMeshSamplingDriver::buildDomainCdsa() {
  const int numdof = domain_->numdof();
  SimpleBuffer<int> bc(numdof);
  SimpleBuffer<double> bcx(numdof);

  domain_->make_bc(bc.array(), bcx.array());
  domain_->make_constrainedDSA(bc.array());
}

void
ExplicitMeshSamplingDriver::readPodBasis(const BasisFileId &basisId,
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
ExplicitMeshSamplingDriver::outputMeshFile(const FileNameInfo &fileInfo, const MeshDesc &mesh) {
  std::ofstream meshOut(getMeshFilename(fileInfo).c_str());
  meshOut << mesh;
}

std::string
ExplicitMeshSamplingDriver::getMeshFilename(const FileNameInfo &fileInfo) {
  return fileInfo.prefix() + ".sampledmesh.inc";
}

void
ExplicitMeshSamplingDriver::buildRestrictedProjection(const NodalRestrictionMapping &restriction,
                                              const VecBasis &originPod, VecBasis &result) {
  restrict_basis(restriction, originPod, result);
  QrPseudoInversion()(result);
}

void
ExplicitMeshSamplingDriver::buildCombinedRestrictedProjection(const NodalRestrictionMapping &restriction,
                                                      const VecBasis &originPod, const VecBasis &targetPod,
                                                      VecBasis &result) {
  VecBasis preProjection;
  buildRestrictedProjection(restriction, originPod, preProjection);
  result.dimensionIs(targetPod.vectorCount(), preProjection.vectorSize());
  combine_projections(targetPod, originPod, preProjection, result);
}
  
void
ExplicitMeshSamplingDriver::outputReducedOperators(const FileNameInfo &fileInfo,
                                           const VecNodeDof6Conversion &conversion,
                                           const NodalRestrictionMapping &restriction,
                                           const VecBasis &forceProjection) {
  outputExtendedPod(BasisFileId(fileInfo, BasisId::FORCE, BasisId::GAPPY_POD), conversion, restriction, forceProjection);
}

void
ExplicitMeshSamplingDriver::outputExtendedPod(const BasisFileId &fileId,
                                      const VecNodeDof6Conversion &conversion,
                                      const NodalRestrictionMapping &restriction,
                                      const VecBasis &basis) {
  BasisOutputStream out(fileId, conversion);
  writeExtendedVectors(out, basis, restriction);
}

void
ExplicitMeshSamplingDriver::outputRestrictedPod(const BasisFileId &fileId,
                                        const VecNodeDof6Conversion &conversion,
                                        const NodalRestrictionMapping &restriction,
                                        const VecBasis &basis) {
  BasisOutputStream out(fileId, conversion);
  writeRestrictedVectors(out, basis, restriction);
}

void
ExplicitMeshSamplingDriver::printReadPodInfo(BasisId::Type type, int vectorCount) {
  std::string name = toString(type);
  assert(!name.empty());
  name[0] = static_cast<char>(std::toupper(name[0]));
  filePrint(stderr, "%s POD subspace of dimension %d\n", name.c_str(), vectorCount);
}

void
ExplicitMeshSamplingDriver::printGappyInfo(const VecBasis &podBasis, const NodalRestrictionMapping &sampling) {
  filePrint(stderr, "Gappy mesh has %d sample nodes\n", sampling.sampleNodeCount());

  const int sampleLocationCount = sampling.restrictedInfo();
  const int podDimensionMax = podBasis.vectorCount();
  const double aspectRatio = static_cast<double>(sampleLocationCount) / podDimensionMax;
  filePrint(stderr, "Gappy aspect ratio = Sampled dofs / POD dimension = %d / %d = %g\n",
                    sampleLocationCount, podDimensionMax, aspectRatio);
}

GenSparseMatrix<double> *
ExplicitMeshSamplingDriver::massMatrixNew() {
  AllOps<double> allOps;
  allOps.M = domain->constructDBSparseMatrix<double>();
  domain->buildOps(allOps, 0.0, 0.0, 0.0, NULL, NULL, NULL, NULL, false);
 
  return allOps.M;
}

} /* end namespace Rom */

Rom::DriverInterface *explicitMeshSamplingDriverNew(Domain *domain) {
  return new Rom::ExplicitMeshSamplingDriver(domain);
}
