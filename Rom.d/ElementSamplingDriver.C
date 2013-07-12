#include "ElementSamplingDriver.h"

#include "VecBasis.h"
#include "BasisOps.h" 
#include "FileNameInfo.h"
#include "BasisFileStream.h"
#include "VecBasisFile.h"
#include "SimpleBuffer.h"
#include "RenumberingUtils.h"
#include "MeshDesc.h"

#include "VecBasisOps.h"

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

void
outputMeshFile(const FileNameInfo &fileInfo, const MeshDesc &mesh, bool firstTime) {
  const std::ios_base::openmode mode = (firstTime) ? std::ios_base::out : std::ios_base::app;  //out = standard output, app = append
  std::ofstream meshOut(getMeshFilename(fileInfo).c_str(), mode);
  meshOut << mesh;
}

void
outputFullWeights(const FileNameInfo &fileInfo, const Vector &weights, const std::vector<int> &elemIds,
                  bool firstTime = true) {
  assert(weights.size() == elemIds.size());

  const std::string fileName = domain->solInfo().reducedMeshFile;
  const std::ios_base::openmode mode = (firstTime) ? std::ios_base::out : std::ios_base::app;
  std::ofstream weightOut(fileName.c_str(), mode);

  if(firstTime) weightOut << "ATTRIBUTES\n";
  for (int i = 0, iEnd = weights.size(); i != iEnd; ++i) {
    weightOut << elemIds[i] + 1 << " 1 " << "HRC" << " " << weights[i] << "\n";
  }
}


// Member functions
// ================
template<typename MatrixBufferType, typename SizeType>
int
ElementSamplingDriver<MatrixBufferType,SizeType>::elementCount() const {
  return domain_->numElements();
}

template<typename MatrixBufferType, typename SizeType>
int
ElementSamplingDriver<MatrixBufferType,SizeType>::vectorSize() const {
  return domain_->numUncon();
}

template<typename MatrixBufferType, typename SizeType>
ElementSamplingDriver<MatrixBufferType,SizeType>::ElementSamplingDriver(Domain *d) :
  SingleDomainDynamic(d),
  domain_(d),
  corotators_(NULL),
  geomState_(NULL),
  kelArray_(NULL),
  melArray_(NULL),
  veloc_(NULL),
  accel_(NULL)
{}

template<typename MatrixBufferType, typename SizeType>
ElementSamplingDriver<MatrixBufferType,SizeType>::~ElementSamplingDriver() {
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

  if(veloc_) delete veloc_;
  if(accel_) delete accel_;
}

template<typename MatrixBufferType, typename SizeType>
void
ElementSamplingDriver<MatrixBufferType,SizeType>::assembleTrainingData(const VecBasis &displac, std::vector<double>::iterator timeStampFirst, const VecBasis &podBasis,
                                                                       typename MatrixBufferType::iterator elemContributions,
                                                                       Vector &trainingTarget, VecBasis *veloc, VecBasis *accel) {
  const int podVectorCount = podBasis.vectorCount();
  const int snapshotCount = displac.vectorCount();

  // Temporary buffers shared by all iterations
  Vector elemTarget(podVectorCount);
  SimpleBuffer<double> elementForce(domain_->maxNumDOF());

  for (int iElem = 0; iElem != elementCount(); ++iElem) {
    filePrint(stderr,"\r %4.2f%% complete", double(iElem)/double(elementCount())*100.);
    std::vector<double>::iterator timeStampIt = timeStampFirst;
    int *nodes = domain_->getElementSet()[iElem]->nodes();
    for (int iSnap = 0; iSnap != snapshotCount; ++iSnap) {
      geomState_->explicitUpdate(domain_->getNodes(), domain_->getElementSet()[iElem]->numNodes(),
                                 nodes, displac[iSnap]); // just set the state at the nodes of element iElem
      if(veloc) geomState_->setVelocity(domain_->getElementSet()[iElem]->numNodes(), nodes,
                                        (*veloc)[iSnap], 2); // just set the velocity at the nodes of element iElem
      if(accel) geomState_->setAcceleration(domain_->getElementSet()[iElem]->numNodes(), nodes,
                                            (*veloc)[iSnap], 2); // just set the acceleration at the nodes of element iElem
      // Evaluate and store element contribution at training configuration
      domain_->getElemInternalForce(*geomState_, *timeStampIt, geomState_, *(corotators_[iElem]), elementForce.array(), kelArray_[iElem]);
      if(domain_->getElementSet()[iElem]->hasRot()) {
        domain_->transformElemStiffAndForce(*geomState_, elementForce.array(), kelArray_[iElem], iElem, false);
        domain_->getElemFictitiousForce(iElem, *geomState_, elementForce.array(), kelArray_[iElem],
                                        *timeStampIt, geomState_, melArray_[iElem], false);
      }

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
  filePrint(stderr,"\r %4.2f%% complete\n", 100.);
}

template<typename MatrixBufferType, typename SizeType>
void
ElementSamplingDriver<MatrixBufferType,SizeType>::solve() {
  Vector solution;
  computeSolution(solution);
  postProcess(solution);
}

template<typename MatrixBufferType, typename SizeType>
void
ElementSamplingDriver<MatrixBufferType,SizeType>::computeSolution(Vector &solution, bool verboseFlag) {
  preProcess();

  const int podVectorCount = podBasis_.vectorCount();
  const int snapshotCount = displac_.vectorCount();

  // DEBUG: Print info
  if(verboseFlag) {
    std::cout << "podVectorCount = " << podVectorCount << ", "
              << "snapshotCount = " << snapshotCount << ", "
              << "elementCount = " << elementCount() << "\n";
  }

  solver_.problemSizeIs(podVectorCount*snapshotCount, elementCount());

  // Training target is the sum of elementary contributions
  Vector trainingTarget(podVectorCount*snapshotCount, 0.0);

  assembleTrainingData(displac_, timeStamps_.begin(), podBasis_, solver_.matrixBuffer(), trainingTarget,
                       veloc_, accel_);

  double targetMagnitude = norm(trainingTarget);

  // Setup and solve optimization problem
  const double relativeTolerance = domain_->solInfo().tolPodRom;
  solver_.relativeToleranceIs(relativeTolerance);
  copy(trainingTarget, solver_.rhsBuffer());

  solver_.verboseFlagIs(verboseFlag);
  solver_.solve();
  
  if(verboseFlag) {
    std::cout << "Primal solution:";
    for (int elemRank = 0; elemRank != elementCount(); ++elemRank) {
      std::cout << " " << solver_.solutionEntry(elemRank);
    }
    std::cout << "\n";

    std::cout << "Error magnitude / Absolute tolerance = " << solver_.errorMagnitude() << " / " << solver_.relativeTolerance() * targetMagnitude << "\n";
    std::cout << "1-norm of primal solution = " << std::accumulate(solver_.solutionBuffer(), solver_.solutionBuffer() + solver_.unknownCount(), 0.0) << "\n";
  }

  // Read solution
  solution.initialize(elementCount());
  std::copy(solver_.solutionBuffer(), solver_.solutionBuffer() + elementCount(), solution.data());
}

template<typename MatrixBufferType, typename SizeType>
void
ElementSamplingDriver<MatrixBufferType,SizeType>::postProcess(Vector &solution, bool firstTime, bool verboseFlag) {

  const FileNameInfo fileInfo;
  std::set<int> sampleElemRanks;
  {
    for (int iElem = 0; iElem != elementCount(); ++iElem) {
      if (solution[iElem] > 0.0) {
        sampleElemRanks.insert(sampleElemRanks.end(), iElem);
      }
    }
  }

  //Element numbering: Packed to input
  std::vector<int> packedToInput(elementCount());
  Elemset &inputElemSet = *(geoSource->getElemSet());
  for (int iElem = 0, iElemEnd = inputElemSet.size(); iElem != iElemEnd; ++iElem) {
    Element *elem = inputElemSet[iElem];
    if (elem) {
      //PJSA const int iPackElem = geoSource->glToPackElem(iElem);
      const int iPackElem = domain_->glToPackElem(iElem);
      assert(iPackElem < packedToInput.size());
      if(iPackElem >= 0) packedToInput[iPackElem] = iElem;
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
  
  const MeshRenumbering meshRenumbering(sampleElemIds.begin(), sampleElemIds.end(), *elemToNode, verboseFlag);
  const MeshDesc reducedMesh(domain_, geoSource, meshRenumbering, weights);
  try {
    outputMeshFile(fileInfo, reducedMesh, firstTime);
  }
  catch(std::exception& e) {
    std::cerr << "caught exception: " << e.what() << endl;
  }
  outputFullWeights(fileInfo, solution, packedToInput, firstTime);
}

template<typename MatrixBufferType, typename SizeType>
void
ElementSamplingDriver<MatrixBufferType,SizeType>::preProcess() {

  domain_->preProcessing();
  buildDomainCdsa();
  domain_->makeAllDOFs();
  
  StaticTimers dummyTimes;
  GenFullSquareMatrix<double> *dummyGeomKelArray = NULL;
  const bool buildMelArray = true;
  domain_->computeGeometricPreStress(corotators_, geomState_, kelArray_, &dummyTimes, dummyGeomKelArray, melArray_, buildMelArray);
  if(domain_->nDirichlet() > 0) {
    geomState_->updatePrescribedDisplacement(domain_->getDBC(), domain_->nDirichlet(), domain_->getNodes());
  }

  const FileNameInfo fileInfo;
  
  // Read order reduction data
  const VecNodeDof6Conversion vecDofConversion(*domain_->getCDSA());
  assert(vectorSize() == vecDofConversion.vectorSize());

  {
    BasisInputStream in(BasisFileId(fileInfo, BasisId::STATE, BasisId::POD), vecDofConversion);
    const int podSizeMax = domain_->solInfo().maxSizePodRom;
    if (podSizeMax != 0) {
      readVectors(in, podBasis_, podSizeMax);
    } else {
      readVectors(in, podBasis_);
    }
  }


// Read state snapshots
  VecBasis snapshots;
  {
    BasisInputStream in(BasisFileId(fileInfo, BasisId::STATE, BasisId::SNAPSHOTS), vecDofConversion);
    const int skipFactor = domain_->solInfo().skipPodRom;
    const int skipOffSet = domain_->solInfo().skipOffSet;
    const int basisStateCount = (in.size() % 2) + (in.size() - skipOffSet) / skipFactor;

    snapshots.dimensionIs(basisStateCount, in.vectorSize());
    timeStamps_.reserve(basisStateCount);

    int count = 0;
    int skipCounter = skipFactor-skipOffSet;
    while (count < basisStateCount) {
      std::pair<double, double *> data;
      data.second = snapshots[count].data();
      in >> data;
      assert(in);
      if (skipCounter == skipFactor) {
        timeStamps_.push_back(data.first);
        skipCounter = 1;
        ++count;
      } else {
        ++skipCounter;
      }
    }

    assert(timeStamps_.size() == basisStateCount);
  }

  // Read velocity snapshots
  VecBasis *velocSnapshots = 0;
  if(domain_->solInfo().velocPodRomFile != "") {
    //std::cerr << "reading velocity snapshots from file " << domain->solInfo().velocPodRomFile << std::endl;
    std::vector<double> timeStamps;
    velocSnapshots = new VecBasis;
    BasisInputStream in(BasisFileId(fileInfo, BasisId::VELOCITY, BasisId::SNAPSHOTS), vecDofConversion);
    const int skipFactor = domain_->solInfo().skipPodRom;
    const int skipOffSet = domain_->solInfo().skipOffSet;
    const int basisStateCount = 1 + (in.size() - 1) / skipFactor;

    velocSnapshots->dimensionIs(basisStateCount, in.vectorSize());
    timeStamps.reserve(basisStateCount);

    int count = 0;
    int skipCounter = skipFactor-skipOffSet;
    while (count < basisStateCount) {
      std::pair<double, double *> data;
      data.second = (*velocSnapshots)[count].data();
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
    // TODO: check that timeStamps for velocity snapshots match state snapshots
  }

  // Read acceleration snapshots
  VecBasis *accelSnapshots = 0;
  if(domain_->solInfo().accelPodRomFile != "") {
    //std::cerr << "reading acceleration snapshots from file " << domain->solInfo().accelPodRomFile << std::endl;
    std::vector<double> timeStamps;
    accelSnapshots = new VecBasis;
    BasisInputStream in(BasisFileId(fileInfo, BasisId::ACCELERATION, BasisId::SNAPSHOTS), vecDofConversion);
    const int skipFactor = domain_->solInfo().skipPodRom;
    const int skipOffSet = domain_->solInfo().skipOffSet;
    const int basisStateCount = 1 + (in.size() - 1) / skipFactor;

    accelSnapshots->dimensionIs(basisStateCount, in.vectorSize());
    timeStamps.reserve(basisStateCount);

    int count = 0;
    int skipCounter = skipFactor-skipOffSet;
    while (count < basisStateCount) {
      std::pair<double, double *> data;
      data.second = (*accelSnapshots)[count].data();
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
    // TODO: check that timeStamps for acceleration snapshots match state snapshots
  }

  const int podVectorCount = podBasis_.vectorCount();
  const int snapshotCount = snapshots.vectorCount();

  // Temporary buffers shared by all iterations
  Vector podComponents(podVectorCount);

  // Project snapshots on POD basis to get training configurations
  displac_.dimensionIs(snapshotCount, vectorSize());
  for (int iSnap = 0; iSnap != snapshotCount; ++iSnap) {
    expand(podBasis_, reduce(podBasis_, snapshots[iSnap], podComponents), displac_[iSnap]);
  }

  if(velocSnapshots) {
    veloc_ = new VecBasis(velocSnapshots->vectorCount(), vectorSize());

    // Project velocity snapshots on POD basis to get training configurations
    for (int iSnap = 0; iSnap != velocSnapshots->vectorCount(); ++iSnap) {
      expand(podBasis_, reduce(podBasis_, (*velocSnapshots)[iSnap], podComponents), (*veloc_)[iSnap]);
    }
    delete velocSnapshots;
  }

  if(accelSnapshots) {
    accel_ = new VecBasis(accelSnapshots->vectorCount(), vectorSize());

    // Project acceleration snapshots on POD basis to get training configurations
    for (int iSnap = 0; iSnap != accelSnapshots->vectorCount(); ++iSnap) {
      expand(podBasis_, reduce(podBasis_, (*accelSnapshots)[iSnap], podComponents), (*accel_)[iSnap]);
    }
    delete accelSnapshots;
  }

  double beta = domain_->solInfo().newmarkBeta;
  if(beta == 0.0) {
     filePrint(stderr,"... Renormalizing Projection Basis ...\n");
     VecBasis normalizedBasis;
     DynamMat * dummyDynOps = SingleDomainDynamic::buildOps(1.0,0.0,0.0);

     assert(dummyDynOps->M);
     const GenSparseMatrix<double> &fullMass = *(dummyDynOps->M);
     renormalized_basis(fullMass, podBasis_, normalizedBasis);
     podBasis_.swap(normalizedBasis);
  }

}

template<typename MatrixBufferType, typename SizeType>
void
ElementSamplingDriver<MatrixBufferType,SizeType>::buildDomainCdsa() {
  const int numdof = domain_->numdof();
  SimpleBuffer<int> bc(numdof);
  SimpleBuffer<double> bcx(numdof);

  domain_->make_bc(bc.array(), bcx.array());
  domain_->make_constrainedDSA(bc.array());
}

} // end namespace Rom

Rom::DriverInterface *elementSamplingDriverNew(Domain *d) {
#ifdef USE_STXXL
  if(d->solInfo().oocPodRom)
    // external vector of double's with 16 blocks per page, the cache with 32 pages, and 8 MB blocks (i.e. total cache is 4GB)
    return new Rom::ElementSamplingDriver<stxxl::VECTOR_GENERATOR<double,16,32,8388608,stxxl::RC,stxxl::random>::result,stxxl::uint64>(d);
  else
#endif
  return new Rom::ElementSamplingDriver<std::vector<double>,size_t>(d);
}
