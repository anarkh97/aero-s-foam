#include "ElementSamplingDriver.h"

#include "VecBasis.h"
#include "BasisOps.h" 
#include "FileNameInfo.h"
#include "BasisFileStream.h"
#include "VecBasisFile.h"
#include "SimpleBuffer.h"
#include "RenumberingUtils.h"
#include "MeshDesc.h"

#include <Driver.d/GeoSource.h>
#include <Driver.d/Domain.h>
#include <Math.d/Vector.h>
#include <Timers.d/StaticTimers.h>
#include <Utils.d/Connectivity.h>
#include <Element.d/Element.h>

#include <cstddef>
#include <algorithm>
#include <numeric>
#include <iterator>
#include <set>
#include <map>
#include <vector>
#include <memory>
#include <fstream>
#include <string>
#include <utility>

#include <cassert>
#include <iostream>

extern GeoSource *geoSource;

namespace Rom {

// Non-member functions
// ====================
template <typename Scalar>
inline
void
copy(const GenVector<Scalar> &v, Scalar *target) {
  const double *originBuf = v.data();
  std::copy(originBuf, originBuf + v.size(), target);
}

template <typename Scalar, typename OutputIterator>
inline
void
copy(const GenVector<Scalar> &v, OutputIterator target) {
  const double *originBuf = v.data();
  std::copy(originBuf, originBuf + v.size(), target);
}

inline
std::string
getMeshFilename(const FileNameInfo &fileInfo) {
  return fileInfo.prefix() + ".elementmesh.inc";
}

inline
void
outputMeshFile(const FileNameInfo &fileInfo, const MeshDesc &mesh) {
  std::ofstream meshOut(getMeshFilename(fileInfo).c_str());
  meshOut << mesh;
}

void
outputFullWeights(const FileNameInfo &fileInfo, const Vector &weights, const std::vector<int> &elemIds) {
  assert(weights.size() == elemIds.size());

  const std::string fileName = domain->solInfo().reducedMeshFile; //fileInfo.prefix() + ".attributes.inc";
  std::ofstream weightOut(fileName.c_str());

  weightOut << "ATTRIBUTES\n";
  for (int i = 0, iEnd = weights.size(); i != iEnd; ++i) {
    weightOut << elemIds[i] + 1 << " 1 " << "HRC" << " " << weights[i] << "\n";
  }
}

// Member functions
// ================
inline
int
ElementSamplingDriver::elementCount() const {
  return domain_->numElements();
}

inline
int
ElementSamplingDriver::vectorSize() const {
  return domain_->numUncon();
}

ElementSamplingDriver::ElementSamplingDriver(Domain *d) :
  domain_(d),
  corotators_(NULL),
  geomState_(NULL),
  kelArray_(NULL)
{}

ElementSamplingDriver::~ElementSamplingDriver() {
  if (corotators_) {
    for (int iElem = 0; iElem != elementCount(); ++iElem) {
      const Corotator * c = corotators_[iElem];
      if (!dynamic_cast<const Element *>(c)) {
        delete corotators_[iElem];
      }
    }
    delete[] corotators_;
  }

  delete geomState_;
}

template <typename DblFwdIt>
void
ElementSamplingDriver::assembleTrainingData(const VecBasis &snapshots, DblFwdIt timeStampFirst, const VecBasis &podBasis,
                                            typename SparseNonNegativeLeastSquaresSolver::MatrixBufferType::iterator elemContributions,
                                            Vector &trainingTarget) {
  const int podVectorCount = podBasis.vectorCount();
  const int snapshotCount = snapshots.vectorCount();

  // Temporary buffers shared by all iterations
  VecBasis displac(snapshotCount, vectorSize());
  Vector podComponents(podVectorCount), elemTarget(podVectorCount);
  SimpleBuffer<double> elementForce(domain_->maxNumDOF());

  // Project snapshots on POD basis to get training configurations
  for (int iSnap = 0; iSnap != snapshotCount; ++iSnap) {
    expand(podBasis, reduce(podBasis, snapshots[iSnap], podComponents), displac[iSnap]);
  }

  for (int iElem = 0; iElem != elementCount(); ++iElem) {
    filePrint(stderr,"\r %4.2f%% complete", double(iElem)/double(elementCount())*100.);
    DblFwdIt timeStampIt = timeStampFirst;
    int *nodes = domain_->getElementSet()[iElem]->nodes();
    for (int iSnap = 0; iSnap != snapshotCount; ++iSnap) {
      //geomState_->explicitUpdate(domain_->getNodes(), displac[iSnap]);
      geomState_->explicitUpdate(domain_->getNodes(), domain_->getElementSet()[iElem]->numNodes(),
                                 nodes, displac[iSnap]); // just update the nodes of element iElem
      // Evaluate and store element contribution at training configuration
      domain_->getElemInternalForce(*geomState_, *timeStampIt, NULL, *(corotators_[iElem]), elementForce.array(), kelArray_[iElem]);
      elemTarget.zero();
      const int dofCount = kelArray_[iElem].dim();
      for (int iDof = 0; iDof != dofCount; ++iDof) {
        const int vecLoc = domain_->getCDSA()->getRCN((*domain_->getAllDOFs())[iElem][iDof]);
        if (vecLoc >= 0) {
          const double dofForce = elementForce[iDof];
          for (int iPod = 0; iPod != podVectorCount; ++iPod) {
            const double contrib = dofForce * podBasis[iPod][vecLoc];
            elemTarget[iPod] += contrib;
          }
        }
      }
      for(int iPod = 0; iPod != podVectorCount; ++iPod) {
        *elemContributions = elemTarget[iPod];
        elemContributions++;
        trainingTarget[podVectorCount * iSnap + iPod] += elemTarget[iPod];
      }
    }
    delete [] nodes;
    timeStampIt++;
  }
 filePrint(stderr,"\n");
}

void
ElementSamplingDriver::solve() {
  preProcess();

  const FileNameInfo fileInfo;
  
  Vector solution;
  {

    // Read order reduction data
    const VecNodeDof6Conversion vecDofConversion(*domain_->getCDSA());
    assert(vectorSize() == vecDofConversion.vectorSize());

    VecBasis podBasis;
    {
      BasisInputStream in(BasisFileId(fileInfo, BasisId::STATE, BasisId::POD), vecDofConversion);

      const int podSizeMax = domain_->solInfo().maxSizePodRom;
      if (podSizeMax != 0) {
        readVectors(in, podBasis, podSizeMax);
      } else {
        readVectors(in, podBasis);
      }
    }

    VecBasis snapshots;
    std::vector<double> timeStamps;
    {
      BasisInputStream in(BasisFileId(fileInfo, BasisId::STATE, BasisId::SNAPSHOTS), vecDofConversion);

      const int skipFactor = domain->solInfo().skipPodRom;
      const int basisStateCount = 1 + (in.size() - 1) / skipFactor;

      snapshots.dimensionIs(basisStateCount, in.vectorSize());
      timeStamps.reserve(basisStateCount);

      int count = 0;
      int skipCounter = skipFactor;
      while (count < basisStateCount) {
        std::pair<double, double *> data;
        data.second = snapshots[count].data();
        in >> data;
        assert(in);
        if (skipCounter == skipFactor) {
          timeStamps.push_back(data.first);
          skipCounter = 1;
          ++count;
        } else {
          ++skipCounter;
        }
      }

      assert(timeStamps.size() == basisStateCount);
    }

    const int podVectorCount = podBasis.vectorCount();
    const int snapshotCount = snapshots.vectorCount();

    // DEBUG: Print info
    std::cout << "podVectorCount = " << podVectorCount << ", "
              << "snapshotCount = " << snapshotCount << ", "
              << "elementCount = " << elementCount() << "\n";

    solver_.problemSizeIs(podVectorCount*snapshotCount, elementCount());

    // Training target is the sum of elementary contributions
    Vector trainingTarget(podVectorCount*snapshotCount, 0.0);

    assembleTrainingData(snapshots, timeStamps.begin(), podBasis, solver_.matrixBuffer(), trainingTarget);

    double targetMagnitude = norm(trainingTarget);

    // Setup and solve optimization problem
    const double relativeTolerance = domain_->solInfo().tolPodRom;
    solver_.relativeToleranceIs(relativeTolerance);
    copy(trainingTarget, solver_.rhsBuffer());

    solver_.solve();
    
    std::cout << "Primal solution:";
    for (int elemRank = 0; elemRank != elementCount(); ++elemRank) {
      std::cout << " " << solver_.solutionEntry(elemRank);
    }
    std::cout << "\n";

    std::cout << "Error magnitude / Absolute tolerance = " << solver_.errorMagnitude() << " / " << solver_.relativeTolerance() * targetMagnitude << "\n";
    std::cout << "1-norm of primal solution = " << std::accumulate(solver_.solutionBuffer(), solver_.solutionBuffer() + solver_.unknownCount(), 0.0) << "\n";

    // Read solution
    solution.initialize(elementCount());
    std::copy(solver_.solutionBuffer(), solver_.solutionBuffer() + elementCount(), solution.data());
  }

  std::set<int> sampleElemRanks;
  {
    for (int iElem = 0; iElem != elementCount(); ++iElem) {
      if (solution[iElem] > 0.0) {
        sampleElemRanks.insert(sampleElemRanks.end(), iElem);
      }
    }
  }

  // Element numbering: Packed to input
  std::vector<int> packedToInput(elementCount());
  Elemset &inputElemSet = *(geoSource->getElemSet());
  for (int iElem = 0, iElemEnd = inputElemSet.size(); iElem != iElemEnd; ++iElem) {
    Element *elem = inputElemSet[iElem];
    if (elem) {
      const int iPackElem = geoSource->glToPackElem(iElem);
      assert(iPackElem >= 0 && iPackElem < packedToInput.size());
      packedToInput[iPackElem] = iElem;
    }
  }

  std::vector<int> sampleElemIds;
  sampleElemIds.reserve(sampleElemRanks.size());
  std::map<int, double> weights;
  for (std::set<int>::const_iterator it = sampleElemRanks.begin(), it_end = sampleElemRanks.end(); it != it_end; ++it) {
    const int elemRank = packedToInput[*it];
    weights.insert(std::make_pair(elemRank, solution[*it]));
    sampleElemIds.push_back(elemRank);
  }

  std::auto_ptr<Connectivity> elemToNode(new Connectivity(&inputElemSet));
  const MeshRenumbering meshRenumbering(sampleElemIds.begin(), sampleElemIds.end(), *elemToNode);

  const MeshDesc reducedMesh(domain_, geoSource, meshRenumbering, weights);
  outputMeshFile(fileInfo, reducedMesh);

  outputFullWeights(fileInfo, solution, packedToInput);
}

void
ElementSamplingDriver::preProcess() {
  domain_->preProcessing();
  buildDomainCdsa();
  domain_->makeAllDOFs();
  
  StaticTimers dummyTimes;
  GenFullSquareMatrix<double> *dummyGeomKelArray = NULL;
  GenFullSquareMatrix<double> *dummyMelArray = NULL;
  const bool buildMelArray = false;
  domain_->computeGeometricPreStress(corotators_, geomState_, kelArray_, &dummyTimes, dummyGeomKelArray, dummyMelArray, buildMelArray);
  if(domain_->nDirichlet() > 0) {
    geomState_->updatePrescribedDisplacement(domain_->getDBC(), domain_->nDirichlet(), domain_->getNodes());
  }
}

void
ElementSamplingDriver::buildDomainCdsa() {
  const int numdof = domain_->numdof();
  SimpleBuffer<int> bc(numdof);
  SimpleBuffer<double> bcx(numdof);

  domain_->make_bc(bc.array(), bcx.array());
  domain_->make_constrainedDSA(bc.array());
}

} // end namespace Rom

Rom::DriverInterface *elementSamplingDriverNew(Domain *d) {
  return new Rom::ElementSamplingDriver(d);
}
