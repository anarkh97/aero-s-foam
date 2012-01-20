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

  const std::string fileName = fileInfo.prefix() + ".fullweights.inc";
  std::ofstream weightOut(fileName.c_str());

  weightOut << "ELLUMP\n";
  for (int i = 0, iEnd = weights.size(); i != iEnd; ++i) {
    weightOut << elemIds[i] + 1 << " " << weights[i] << "\n";
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
ElementSamplingDriver::assembleTrainingData(const VecBasis &snapshots, DblFwdIt timeStampFirst, const VecBasis &podBasis, VecBasis &elemContributions) {
  const int podVectorCount = podBasis.vectorCount();
  const int snapshotCount = snapshots.vectorCount();

  elemContributions.dimensionIs(elementCount(), snapshotCount * podVectorCount);
  std::fill_n(&elemContributions[0][0], elemContributions.vectorSize() * elemContributions.vectorCount(), 0.0);

  // Temporary buffers shared by all iterations
  Vector displac(vectorSize());
  Vector podComponents(podVectorCount);
  SimpleBuffer<double> elementForce(domain_->maxNumDOF());

  DblFwdIt timeStampIt = timeStampFirst;
  for (int iSnap = 0; iSnap != snapshotCount; ++iSnap) {
    // Project snapshot on POD basis to get training configuration
    expand(podBasis, reduce(podBasis, snapshots[iSnap], podComponents), displac);
    geomState_->explicitUpdate(domain_->getNodes(), displac);
    // Evaluate and store elementary contributions at training configuration
    for (int iElem = 0; iElem != elementCount(); ++iElem) {
      domain_->getElemStiffAndForce(*geomState_, *timeStampIt, NULL, *(corotators_[iElem]), elementForce.array(), kelArray_[iElem]);
      double * const targetBuffer = elemContributions[iElem].data() + (podVectorCount * iSnap);
      const int dofCount = kelArray_[iElem].dim();
      for (int iDof = 0; iDof != dofCount; ++iDof) {
        const int vecLoc = domain_->getCDSA()->getRCN((*domain_->getAllDOFs())[iElem][iDof]);
        if (vecLoc >= 0) {
          const double dofForce = elementForce[iDof];
          for (int iPod = 0; iPod != podVectorCount; ++iPod) {
            const double contrib = dofForce * podBasis[iPod][vecLoc];
            targetBuffer[iPod] += contrib;
          }
        }
      }
    }
    timeStampIt++;
  }
}

void
ElementSamplingDriver::solve() {
  preProcess();

  const FileNameInfo fileInfo;
  
  Vector solution;
  {
    // Training data
    VecBasis elemContributions;
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

      assembleTrainingData(snapshots, timeStamps.begin(), podBasis, elemContributions);
    }

    double targetMagnitude;
    {
      // Training target is the sum of elementary contributions
      Vector trainingTarget(elemContributions.vectorSize(), 0.0);
      for (VecBasis::const_iterator it = elemContributions.begin(), it_end = elemContributions.end(); it != it_end; ++it) {
        trainingTarget += *it;
      }
      targetMagnitude = norm(trainingTarget);

      // Setup and solve optimization problem
      const double relativeTolerance = domain_->solInfo().tolPodRom;
      solver_.relativeToleranceIs(relativeTolerance);
      solver_.problemSizeIs(elemContributions.vectorSize(), elemContributions.vectorCount());
      {
        int col = 0;
        for (VecBasis::const_iterator it = elemContributions.begin(); it != elemContributions.end(); ++it) {
          copy(*it, solver_.matrixColBuffer(col));
          ++col;
        }
      }
      copy(trainingTarget, solver_.rhsBuffer());
    }

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
    geomState_->updatePrescribedDisplacement(domain_->getDBC(), domain_->nDirichlet());
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
