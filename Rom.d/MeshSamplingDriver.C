#include "MeshSamplingDriver.h"

#include "GreedyUtils.h"
#include "ConnectivityUtils.h"
#include "RenumberingUtils.h"
#include "MeshUtils.h"
#include "VecBasisUtils.h"
#include "VecBasisFile.h"
#include "FileNameInfo.h"

#include "SamplingErrorEvaluation.h"

#include "BasisFileStream.h"
#include "VecNodeDof6Conversion.h"
#include "NodalRestrictionMapping.h"
#include "QrPseudoInversion.h"

#include "MeshOutput.h"

#include "SimpleBuffer.h"

#include <Driver.d/Domain.h>
#include <Driver.d/GeoSource.h> 
#include <Utils.d/DistHelper.h>

#include <map>
#include <vector>
#include <iterator>
#include <cstddef>
#include <fstream>
#include <cmath>

// DEBUG
#include <iostream>

extern GeoSource *geoSource;

MeshSamplingDriver::MeshSamplingDriver(Domain *domain) :
  domain_(domain)
{}

void
MeshSamplingDriver::solve() {
  preProcess();
  
  std::vector<int> sampleNodes = determineSampleNodes();
  buildReducedMesh(sampleNodes);
}

std::vector<int>
MeshSamplingDriver::determineSampleNodes() {
  VecNodeDof6Conversion vecDofConversion(*domain->getCDSA());
  FileNameInfo fileInfo;
  
  // Read bases from file
  std::vector<VecBasis> podBasis(2);
  BasisId::Type types[2] = { BasisId::RESIDUAL, BasisId::JACOBIAN };

  for (BasisId::Type *it = types; it != types + 2; ++it) {
    BasisInputStream input(fileInfo.fileName(BasisId(*it, BasisId::POD)), vecDofConversion);

    const int projectionSubspaceSize = domain->solInfo().maxSizePodRom ?
                                       std::min(domain->solInfo().maxSizePodRom, input.size()) :
                                       input.size();

    readVectors(input, podBasis[it - types], projectionSubspaceSize);
  }

  // Greedy iterations parameters
  const double aspectRatio = domain->solInfo().aspectRatioPodRom;
  const int podVectorCountMax = std::max(podBasis[0].numVec(), podBasis[1].numVec());
  const int targetSampleSize = std::max(static_cast<int>(aspectRatio * podVectorCountMax), podVectorCountMax);
  const int dofsPerNodeEstimate = 6;

  std::vector<int> sampleNodes;
  std::vector<int> sampleLocations;
  
  std::vector<Vector> reconstructionError(2);
  SamplingErrorEvaluation errorEval;

  // Start greedy iterations
  int handledVectorCount = 0;
  int greedyIter = 0;
  const int podVectorCountMin = std::min(podBasis[0].numVec(), podBasis[1].numVec());
  while (sampleLocations.size() < targetSampleSize && handledVectorCount < podVectorCountMin) {
    std::cerr << sampleLocations.size() << " sample locations at greedy iteration: " << greedyIter++ << std::endl;

    // Update the reconstruction error
    handledVectorCount++;
    std::vector<Vector>::iterator errorIt = reconstructionError.begin();
    for (std::vector<VecBasis>::const_iterator podIt = podBasis.begin(); podIt != podBasis.end(); ++podIt) {
      errorEval(podIt->begin(), podIt->begin() + handledVectorCount,
                sampleLocations.begin(), sampleLocations.end(),
                *errorIt++);
    }

    // Zero out the entries corresponding to the already selected indices to avoid selecting them again
    for (std::vector<Vector>::iterator vecIt = reconstructionError.begin(); vecIt != reconstructionError.end(); ++vecIt) {
      Vector &vec = *vecIt;
      for (std::vector<int>::const_iterator it = sampleLocations.begin(); it != sampleLocations.end(); ++it) {
        vec[*it] = 0.0;
      }
    }

    // Greedily determine the best entries and find the corresponding nodes
    const int remainingLocations = targetSampleSize - sampleLocations.size();
    const int remainingPodVectors = (podVectorCountMin - handledVectorCount) + 1;
    const double ratio = static_cast<double>(remainingLocations) / (remainingPodVectors * dofsPerNodeEstimate);
    const int locationsToSelect = static_cast<int>(std::ceil(ratio));

    std::vector<int> selectedLocations;
    max_magnitude_indices(reconstructionError.begin(), reconstructionError.end(),
                          std::back_inserter(selectedLocations), locationsToSelect);

    std::set<int> selectedNodes; // To ensure unicity, since several locations can correspond to one node
    for (std::vector<int>::const_iterator it = selectedLocations.begin(); it != selectedLocations.end(); ++it) {
      selectedNodes.insert(vecDofConversion.nodeDof(*it).nodeRank);
    }

    // Add all the dofs attached to the newly selected nodes
    std::cerr << "Selected nodes: ";
    for (std::set<int>::const_iterator it = selectedNodes.begin(); it != selectedNodes.end(); ++it) {
      std::cerr << *it + 1 << " ";
      sampleNodes.push_back(*it); // Newly selected nodes are not already selected
      vecDofConversion.locations(*it, std::back_inserter(sampleLocations));
    }
    std::cerr << std::endl;
  }

  filePrint(stderr, "Selected %d sample nodes corresponding to %d degrees of freedom\n", sampleNodes.size(), sampleLocations.size());
  filePrint(stderr, "Actual gappy POD aspect ratio = %g\n", static_cast<double>(sampleLocations.size()) / podVectorCountMax);

  // Output all sample nodes according to the original numbering
  // TODO: Consider removal
  const std::string sampleOutName = fileInfo.prefix() + ".samplenodes.check";
  std::ofstream sampleOut(sampleOutName.c_str());
  sampleOut << make_section(sampleNodes.begin(), sampleNodes.end(), SampleNodeTag());

  return sampleNodes;
}

void
MeshSamplingDriver::buildReducedMesh(const std::vector<int> &sampleNodes) {
  std::cerr << "Sample nodes = ";
  for (std::vector<int>::const_iterator it = sampleNodes.begin(); it != sampleNodes.end(); ++it) {
    std::cerr << *it + 1 << " ";
  }
  
  std::cerr << std::endl;
  std::vector<int> reducedNodes;
  neighborhood(*domain->getNodeToNode(), sampleNodes.begin(), sampleNodes.end(), std::back_inserter(reducedNodes));

  std::cerr << "Nodes = ";
  for (std::vector<int>::const_iterator it = reducedNodes.begin(); it != reducedNodes.end(); ++it) {
    std::cerr << *it + 1 << " ";
  }
  std::cerr << std::endl;

  std::map<int, int> nodeRenumbering;
  inverse_numbering(reducedNodes.begin(), reducedNodes.end(), std::inserter(nodeRenumbering, nodeRenumbering.end()));

  std::vector<int> reducedSampleNodes;
  renumber(nodeRenumbering, sampleNodes.begin(), sampleNodes.end(), std::back_inserter(reducedSampleNodes));

  std::cerr << "Node renumbering = ";
  for (std::map<int, int>::const_iterator it = nodeRenumbering.begin(); it != nodeRenumbering.end(); ++it) {
    std::cerr << it->first + 1 << "->" << it->second + 1 << " ";
  }
  std::cerr << std::endl;

  std::vector<int> reducedElems;
  connections(*domain->getNodeToElem(), sampleNodes.begin(), sampleNodes.end(), std::back_inserter(reducedElems));

  std::cerr << "Elements = ";
  for (std::vector<int>::const_iterator it = reducedElems.begin(); it != reducedElems.end(); ++it) {
    std::cerr << *it + 1 << " ";
  }
  std::cerr << std::endl;

  std::map<int, int> elemRenumbering;
  inverse_numbering(reducedElems.begin(), reducedElems.end(), std::inserter(elemRenumbering, elemRenumbering.end()));
  
  std::cerr << "Element renumbering = ";
  for (std::map<int, int>::const_iterator it = elemRenumbering.begin(); it != elemRenumbering.end(); ++it) {
    std::cerr << it->first + 1 << "->" << it->second + 1 << " ";
  }
  std::cerr << std::endl;

  // CoordSet
  CoordSet reducedNodeSet;
  reduce(reducedNodes.begin(), reducedNodes.end(), domain->getNodes(), reducedNodeSet);

  // Elemset
  Elemset reducedElemSet;
  reduce(reducedElems.begin(), reducedElems.end(), nodeRenumbering, *domain->getEset(), reducedElemSet);

  // Attributes
  std::vector<Attrib> reducedAttributes;
  const std::map<int, Attrib> &attribMap = geoSource->getAttributes();
  reduce(elemRenumbering, make_value_iterator(attribMap.begin()), make_value_iterator(attribMap.end()), std::back_inserter(reducedAttributes)); 

  // Boundary conditions
  std::vector<BCond> reducedDirichletBc;
  reduce(nodeRenumbering, domain->getDBC(), domain->getDBC() + domain->nDirichlet(), std::back_inserter(reducedDirichletBc));
  std::vector<BCond> reducedNeumannBc;
  reduce(nodeRenumbering, domain->getNBC(), domain->getNBC() + domain->nNeumann(), std::back_inserter(reducedNeumannBc));

  // Initial conditions
  std::vector<BCond> reducedInitDisp;
  reduce(nodeRenumbering, domain->getInitDisp(), domain->getInitDisp() + domain->numInitDisp(), std::back_inserter(reducedInitDisp));
  reduce(nodeRenumbering, domain->getInitDisp6(), domain->getInitDisp6() + domain->numInitDisp6(), std::back_inserter(reducedInitDisp));

  std::vector<BCond> reducedInitVel;
  reduce(nodeRenumbering, domain->getInitVelocity(), domain->getInitVelocity() + domain->numInitVelocity(), std::back_inserter(reducedInitVel));

  // EFrames
  std::vector<EFrameData> reducedEFrames;
  copy_eframes(reducedElemSet, std::back_inserter(reducedEFrames));

  // Output mesh
  FileNameInfo fileInfo;
  const std::string meshOutName = fileInfo.prefix() + ".sampledmesh.inc";
  std::ofstream meshOut(meshOutName.c_str());

  meshOut << reducedNodeSet;
  meshOut << reducedElemSet;
  if (!reducedEFrames.empty()) {
    meshOut << make_section(reducedEFrames.begin(), reducedEFrames.end());
  }

  meshOut << make_section(reducedAttributes.begin(), reducedAttributes.end());
  meshOut << geoSource->getStructProps();
  
  meshOut << make_section(reducedDirichletBc.begin(), reducedDirichletBc.end(), BCond::Displacements);
  meshOut << make_section(reducedNeumannBc.begin(), reducedNeumannBc.end(), BCond::Forces);
  
  meshOut << make_section(reducedInitDisp.begin(), reducedInitDisp.end(), BCond::Idisplacements);
  meshOut << make_section(reducedInitVel.begin(), reducedInitVel.end(), BCond::Ivelocities);
  
  meshOut << make_section(reducedSampleNodes.begin(), reducedSampleNodes.end(), SampleNodeTag());

  // State POD on reduced mesh
  VecNodeDof6Conversion vecDofConversion(*domain->getCDSA());
  NodalRestrictionMapping reducedMeshMapping(*domain->getCDSA(), reducedNodes.begin(), reducedNodes.end());
  
  BasisInputStream stateInput(fileInfo.fileName(BasisId(BasisId::STATE, BasisId::POD)), vecDofConversion);
  
  const int projectionSubspaceSize = domain->solInfo().maxSizePodRom ?
                                     std::min(domain->solInfo().maxSizePodRom, stateInput.size()) :
                                     stateInput.size();
  filePrint(stderr, "State POD subspace dimension = %d\n", projectionSubspaceSize);
  
  // TODO: factor out 
  VecBasis statePod(projectionSubspaceSize, stateInput.vectorSize());
  readVectors(stateInput, statePod.begin(), statePod.end());

  VecBasis reducedStatePod(statePod.vectorCount(), reducedMeshMapping.restrictedInfo());
  {
    VecBasis::iterator jt = reducedStatePod.begin();
    for (VecBasis::const_iterator it = statePod.begin(); it != statePod.end(); ++it) {
      reducedMeshMapping.restriction(*it, *jt++);
    }
  }
  
  // Reduced DofSetArray
  DofSetArray reducedDsa(reducedNodes.size(), reducedElemSet);
  ConstrainedDSA reducedCdsa(reducedDsa, reducedDirichletBc.size(), reducedDirichletBc.empty() ? NULL : &reducedDirichletBc[0]);
  VecNodeDof6Conversion reducedVecDofConversion(reducedCdsa);

  // Output state POD on reduced mesh
  BasisOutputStream reducedStateOutput(fileInfo.fileName(BasisId(BasisId::STATE, BasisId::GAPPY_POD)), reducedVecDofConversion);
  reducedStateOutput << reducedStatePod;
 
  // Sampled operators A and B 
  NodalRestrictionMapping sampleMapping(*domain->getCDSA(),
                                        sampleNodes.begin(),
                                        sampleNodes.end());

  // Jacobian -> A
  BasisInputStream jacobianInput(fileInfo.fileName(BasisId(BasisId::JACOBIAN, BasisId::POD)), vecDofConversion);

  const int jacobianSubspaceSize = domain->solInfo().maxSizePodRom ?
                                   std::min(domain->solInfo().maxSizePodRom, jacobianInput.size()) :
                                   jacobianInput.size();
  filePrint(stderr, "Jacobian POD subspace dimension = %d\n", jacobianSubspaceSize);
 
  // TODO: factor out 
  VecBasis jacobianPod(jacobianSubspaceSize, stateInput.vectorSize());
  readVectors(jacobianInput, jacobianPod.begin(), jacobianPod.end());
  
  VecBasis jacobianProjection(jacobianPod.vectorCount(), sampleMapping.restrictedInfo());
  {
    VecBasis::iterator jt = jacobianProjection.begin();
    for (VecBasis::const_iterator it = jacobianPod.begin(); it != jacobianPod.end(); ++it) {
      sampleMapping.restriction(*it, *jt++);
    }
  }

  QrPseudoInversion()(jacobianProjection);
 
  NodalRestrictionMapping reducedSampleMapping(reducedCdsa, reducedSampleNodes.begin(), reducedSampleNodes.end());
  Vector buffer(reducedSampleMapping.originInfo());

  BasisOutputStream jacobianProjOutput(fileInfo.fileName(BasisId(BasisId::JACOBIAN, BasisId::GAPPY_POD)), reducedVecDofConversion);
  for (VecBasis::const_iterator vecIt = jacobianProjection.begin(); vecIt != jacobianProjection.end(); ++vecIt) {
    reducedSampleMapping.extension(*vecIt, buffer);
    jacobianProjOutput << buffer;
  }
  
  // Residual, Jacobian -> B
  BasisInputStream residualInput(fileInfo.fileName(BasisId(BasisId::RESIDUAL, BasisId::POD)), vecDofConversion);
  const int residualSubspaceSize = domain->solInfo().maxSizePodRom ?
                                   std::min(domain->solInfo().maxSizePodRom, residualInput.size()) :
                                   residualInput.size();

  filePrint(stderr, "Residual POD subspace dimension = %d\n", residualSubspaceSize);
  
  // TODO: factor out 
  VecBasis residualPod(residualSubspaceSize, stateInput.vectorSize());
  readVectors(residualInput, residualPod.begin(), residualPod.end());
  
  VecBasis residualPreProjection(residualPod.vectorCount(), sampleMapping.restrictedInfo());
  {
    VecBasis::iterator jt = residualPreProjection.begin();
    for (VecBasis::const_iterator it = residualPod.begin(); it != residualPod.end(); ++it) {
      sampleMapping.restriction(*it, *jt++);
    }
  }

  QrPseudoInversion()(residualPreProjection);

  VecBasis residualProjection(jacobianPod.vectorCount(), sampleMapping.restrictedInfo());
  combine_projections(jacobianPod, residualPod, residualPreProjection, residualProjection);
  
  BasisOutputStream residualProjOutput(fileInfo.fileName(BasisId(BasisId::RESIDUAL, BasisId::GAPPY_POD)), reducedVecDofConversion);
  for (VecBasis::const_iterator vecIt = residualProjection.begin(); vecIt != residualProjection.end(); ++vecIt) {
    reducedSampleMapping.extension(*vecIt, buffer);
    residualProjOutput << buffer;
  }
  
}

void
MeshSamplingDriver::preProcess() {
  domain_->preProcessing();
 
  // Build the constrained DofSetArray incorporating the boundary conditions 
  const int numdof = domain_->numdof();
  SimpleBuffer<int> bc(numdof);
  SimpleBuffer<double> bcx(numdof);

  domain_->make_bc(bc.array(), bcx.array());
  domain_->make_constrainedDSA(bc.array());
}

