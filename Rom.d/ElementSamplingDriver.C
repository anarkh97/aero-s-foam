#include "ElementSamplingDriver.h"

#include "DistrBasisFile.h"
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
#include <Driver.d/SysState.h>
#include <Math.d/Vector.h>
#include <Math.d/DBSparseMatrix.h>
#include <Math.d/EiSparseMatrix.h>
#include <Math.d/DiagMatrix.h>
#include <Timers.d/StaticTimers.h>
#include <Utils.d/Connectivity.h>
#include <Element.d/Element.h>
#include <Utils.d/Conwep.d/BlastLoading.h>
#include <Corotational.d/Corotator.h>

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

std::string
getMeshFilename(const FileNameInfo &fileInfo) {
//  return fileInfo.prefix() + ".elementmesh.inc";
  std::string FName(domain->solInfo().reducedMeshFile, std::strlen(domain->solInfo().reducedMeshFile));
  return FName + ".elementmesh.inc";
}

void
outputMeshFile(const FileNameInfo &fileInfo, const MeshDesc &mesh, const int podVectorCount) {
  const std::ios_base::openmode mode = std::ios_base::out; 
  std::ofstream meshOut(getMeshFilename(fileInfo).c_str(), mode);
  filePrint(stderr," ... Writing Mesh File to %s ...\n", getMeshFilename(fileInfo).c_str());
  meshOut.precision(std::numeric_limits<double>::digits10+1);
  std::string basisfile = getMeshFilename(fileInfo).c_str();
  basisfile.append(".compressed.basis");
  if(domain->numInitDisp() || domain->numInitDisp6() || domain->numInitVelocity()) {
    std::string basisfile2(basisfile);
    if(domain->solInfo().useMassNormalizedBasis || domain->solInfo().newmarkBeta == 0) basisfile2.append(".normalized");
    meshOut << "READMODE \"" << basisfile << "\" \"" << basisfile2 << "\" " << podVectorCount << "\n";
  }
  else {
    meshOut << "READMODE \"" << basisfile << "\" " << podVectorCount << "\n";
  }
  if(!domain->solInfo().useMassNormalizedBasis && domain->solInfo().newmarkBeta != 0)
    meshOut << "use_mass_normalized_basis off\n";
  meshOut << "*\n";
  meshOut << mesh;
}

void
outputMeshFile(const FileNameInfo &fileInfo, const MeshDesc &mesh, const std::vector<int> &localBasisSize) {
  const std::ios_base::openmode mode = std::ios_base::out;
  std::ofstream meshOut(getMeshFilename(fileInfo).c_str(), mode);
  filePrint(stderr," ... Writing Mesh File to %s ...\n", getMeshFilename(fileInfo).c_str());
  meshOut.precision(std::numeric_limits<double>::digits10+1);
  std::string basisfile = getMeshFilename(fileInfo).c_str();
  basisfile.append(".compressed.basis");
  if(domain->numInitDisp() || domain->numInitDisp6() || domain->numInitVelocity()) {
    std::cerr << "Error: IDISP/IVEL is not supported yet for local bases\n";
    exit(-1);
  }
  else {
    for(int j=0; j<localBasisSize.size(); ++j) {
      meshOut << "READMODE \"" << basisfile << j+1 << "\" " << localBasisSize[j] << "\n";
    }
  }
  if(!domain->solInfo().useMassNormalizedBasis && domain->solInfo().newmarkBeta != 0)
    meshOut << "use_mass_normalized_basis off\n";
  meshOut << "*\n";
  meshOut << mesh;
}


template<typename WeightsVecType, typename ElemIdsVecType>
void
outputFullWeights(const WeightsVecType &weights, const ElemIdsVecType &elemIds, int j)
{
  assert(weights.size() == elemIds.size());

  std::string fileName = domain->solInfo().reducedMeshFile;
  if(j > -1) {
    std::ostringstream ss;
    ss << ".cluster" << j+1;
    fileName.append(ss.str());
  }
  const std::ios_base::openmode mode = std::ios_base::out;
  std::ofstream weightOut(fileName.c_str(), mode);
  weightOut.precision(std::numeric_limits<double>::digits10+1);
  bool firstTime = true;

  std::map<int, Attrib> &attrib = geoSource->getAttributes();
  int na = geoSource->getNumAttributes();
  int nMaxEle = geoSource->getElemSet()->last();

  std::map<int, Attrib>::iterator *elemAttrib = new std::map<int, Attrib>::iterator[nMaxEle];
  for(int i = 0; i < nMaxEle; ++i) elemAttrib[i] = attrib.end();
  for (std::map<int, Attrib>::iterator it = attrib.begin(); it != attrib.end(); ++it) {
    if(it->second.nele < nMaxEle)
      elemAttrib[it->second.nele] = it;
  }

  weightOut << "ATTRIBUTES\n";
  for (int i = 0, iEnd = weights.size(); i != iEnd; ++i) {
    if(elemAttrib[elemIds[i]] != attrib.end()) {
      // element has an attribute
      Attrib &a = elemAttrib[elemIds[i]]->second;
      weightOut << elemIds[i]+1 << " " << a.attr+1 << " ";
      if(a.cmp_attr >= 0) {
        if(a.cmp_frm >= -0) weightOut << a.cmp_attr+1 << " " << a.cmp_frm+1 << " ";
        else  weightOut << a.cmp_attr+1 << " THETA " << a.cmp_theta << " ";
      }
    }
    else {
      // element has no attribute, however we still need to define the weight for it
      weightOut << elemIds[i]+1 << " ";
    }
    if(domain->solInfo().reduceFollower && firstTime) {
      weightOut << "HRC " << weights[i] << " EXTFOL\n";
      firstTime = false;
    }
    else {
      weightOut << "HRC " << weights[i] << "\n";
    }
  }

  delete [] elemAttrib;
  weightOut.close();
}

int snapSize(BasisId::Type type, std::vector<int> &snapshotCounts, int j=-1)
{
  // compute the number of snapshots that will be used for a given skipping strategy
  FileNameInfo fileInfo;
  snapshotCounts.clear();
  const int skipFactor = std::max(domain->solInfo().skipPodRom, 1); // skipFactor must be >= 1
  const int skipOffSet = std::max(domain->solInfo().skipOffSet, 0); // skipOffSet must be >= 0
  for(int i = 0; i < FileNameInfo::size(type, BasisId::SNAPSHOTS); i++) {
    std::string fileName = BasisFileId(fileInfo, type, BasisId::SNAPSHOTS, i);
    DistrBasisInputFile in(fileName);
    const double N = in.stateCount() - skipOffSet;
    const int singleSnapshotCount = (N > 0) ? 1+(N-1)/skipFactor : 0;
    snapshotCounts.push_back(singleSnapshotCount);
  }
  // return the total
  return (j<0) ? std::accumulate(snapshotCounts.begin(), snapshotCounts.end(), 0) : snapshotCounts[j];
}

void readAndProjectSnapshots(BasisId::Type type, const int vectorSize, VecBasis &podBasis,
                             const VecNodeDof6Conversion &vecDofConversion,
                             std::vector<int> &snapshotCounts, std::vector<double> &timeStamps, VecBasis &config,
                             SparseMatrix *M, int j = -1)
{
  // if j == -1 then read all of the snapshot files
  // if j >= 0  the read only the (j+1)-th snapshot file
#ifdef PRINT_ESTIMERS
  double t1 = getTime();
#endif
  const int snapshotCount = snapSize(type, snapshotCounts, j);
  filePrint(stderr, " ... Reading in and Projecting %d %s Snapshots ...\n", snapshotCount, toString(type).c_str());

  config.dimensionIs(snapshotCount, vectorSize);
  timeStamps.clear();
  timeStamps.reserve(snapshotCount);

  const int skipFactor = std::max(domain->solInfo().skipPodRom, 1); // skipFactor must be >= 1
  const int skipOffSet = std::max(domain->solInfo().skipOffSet, 0); // skipOffSet must be >= 0
  const int podVectorCount = podBasis.vectorCount();
  Vector snapshot(vectorSize), Msnapshot(vectorSize);
  Vector podComponents(podVectorCount);
  const FileNameInfo fileInfo;

  int offset = 0;
  for(int i = ((j<0)?0:j); i < ((j<0)?FileNameInfo::size(type, BasisId::SNAPSHOTS):(j+1)); i++) {
    std::string fileName = BasisFileId(fileInfo, type, BasisId::SNAPSHOTS, i);
    filePrint(stderr, " ... Processing File: %s ...\n", fileName.c_str());
    BasisInputStream<6> in(fileName, vecDofConversion);

    double s0 = -getTime(), s1 = -51, s2 = 0;
    int count = 0;
    int skipCounter = skipFactor - skipOffSet;
    std::pair<double, double *> data;
    while(count < snapshotCounts[i]) {
      if(skipCounter == skipFactor) {
        data.second = snapshot.data();
        in >> data;
        assert(in);
        if(domain->solInfo().useMassOrthogonalProjection) {
          M->mult(snapshot, Msnapshot);
          expand(podBasis, reduce(podBasis, Msnapshot, podComponents), config[offset+count]);
        }
        else {
          expand(podBasis, reduce(podBasis, snapshot, podComponents), config[offset+count]);
        }
        timeStamps.push_back(data.first);
        skipCounter = 1;
        ++count;
        if((s2-s1 > 50)) { // only print to the screen every 50 milliseconds, otherwise it's too slow...
          s1 = s2;
          filePrint(stderr, "\r ... timeStamp = %8.2e, %3d%% done ...", data.first, (count*100)/snapshotCounts[i]);
        }
        s2 = s0+getTime();
      }
      else {
        in.file().currentStateIndexInc();
        ++skipCounter;
      }
    }

    filePrint(stderr, "\r ... timeStamp = %8.2e, %3d%% done... \n", data.first, 100);
    offset += snapshotCounts[i];
  }
  
  assert(timeStamps.size() == snapshotCount);
#ifdef PRINT_ESTIMERS
  fprintf(stderr, "time for readAndProjectSnapshots = %f\n", (getTime()-t1)/1000.0);
#endif
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
  clean();
}

template<typename MatrixBufferType, typename SizeType>
void
ElementSamplingDriver<MatrixBufferType,SizeType>::clean() {
  if (corotators_) {
    for (int iElem = 0; iElem != elementCount(); ++iElem) {
      const Corotator * c = corotators_[iElem];
      if (!dynamic_cast<const Element *>(c)) {
        delete corotators_[iElem];
      }
    }
    delete[] corotators_;
    corotators_ = NULL;
  }

  if(geomState_) { delete geomState_; geomState_ = NULL; }
  if(kelArray_) { delete [] kelArray_; kelArray_ = NULL; }
  if(melArray_) { delete [] melArray_; melArray_ = NULL; }

  if(veloc_) { delete veloc_; veloc_ = NULL; }
  if(accel_) { delete accel_; accel_ = NULL; }

  timeStamps_.clear();
}

template<typename MatrixBufferType, typename SizeType>
template<typename VecBasisType>
void
ElementSamplingDriver<MatrixBufferType,SizeType>
::assembleTrainingData(const VecBasisType &podBasis, const int podVectorCount, const VecBasisType &displac,
                       const VecBasisType *veloc, const VecBasisType *accel, int j)
{
  std::vector<double>::iterator timeStampFirst = timeStamps_.begin();
  typename MatrixBufferType::iterator elemContributions = solver_.matrixBuffer();
  double *trainingTarget = solver_.rhsBuffer();

  // Temporary buffers shared by all iterations
  Vector elemTarget(podVectorCount);
  SimpleBuffer<double> elementForce(domain_->maxNumDOF());

  BlastLoading::BlastData *conwep;
  if(domain_->solInfo().conwepConfigurations.empty()) conwep = (domain_->solInfo().ConwepOnOff) ? &BlastLoading::InputFileData : NULL;
  else if(domain_->solInfo().conwepConfigurations.size() < domain_->solInfo().statePodRomFile.size()) {
    filePrint(stderr, " *** ERROR: Must provide one Conwep configuration per state snapshot file\n");
    exit(-1);
  }
  domain_->makeElementAdjacencyLists();

  for (int iElem = 0; iElem != elementCount(); ++iElem) {
    filePrint(stderr,"\r %4.2f%% complete", double(iElem)/double(elementCount())*100.);
    std::vector<double>::iterator timeStampIt = timeStampFirst;
    int *nodes = domain_->getElementSet()[iElem]->nodes();
    int iSnap = 0;
    for(int i = ((j<0)?0:j); i < ((j<0)?snapshotCounts_.size():(j+1)); i++) {
      if(!domain_->solInfo().conwepConfigurations.empty()) {
        conwep = &domain_->solInfo().conwepConfigurations[i];
      }
      for (int jSnap = 0; jSnap != snapshotCounts_[i]; ++iSnap, ++jSnap) {
        geomState_->explicitUpdate(domain_->getNodes(), domain_->getElementSet()[iElem]->numNodes(),
            nodes, displac[iSnap]); // just set the state at the nodes of element iElem
        if(veloc) geomState_->setVelocity(domain_->getElementSet()[iElem]->numNodes(), nodes,
            (*veloc)[iSnap], 2); // just set the velocity at the nodes of element iElem
        if(accel) geomState_->setAcceleration(domain_->getElementSet()[iElem]->numNodes(), nodes,
            (*accel)[iSnap], 2); // just set the acceleration at the nodes of element iElem
        // Evaluate and store element contribution at training configuration
        if(corotators_[iElem] && (!domain_->solInfo().getNLInfo().linearelastic ||
          (domain_->getElementSet()[iElem]->isConstraintElement() && domain_->solInfo().getNLInfo().linearelastic == 2))) {
          domain_->getElemInternalForce(*geomState_, *timeStampIt, geomState_, *(corotators_[iElem]), elementForce.array(), kelArray_[iElem]);
        }
        else {
          Vector disp(kelArray_[iElem].dim());
          StackVector force(elementForce.array(), kelArray_[iElem].dim());
          domain->getElementDisp(iElem, *geomState_, disp);
          force.zero();
          kelArray_[iElem].multiply(disp, force, 1.0);
        }
        if(domain_->solInfo().reduceFollower)
          domain_->getElemFollowerForce(iElem, *geomState_, elementForce.array(), elementForce.size(), (corotators_[iElem]),
                                        kelArray_[iElem], 1.0, *timeStampIt, false, conwep);
        if(domain_->getElementSet()[iElem]->hasRot() && !domain_->solInfo().getNLInfo().linearelastic) {
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
        timeStampIt++;
      }
      // reset the element internal states for the next group of snapshots, assuming for now that each group of snapshots
      // was generated by an independent simulation rather than a restart. For a restart it would be approprate to not
      // reset the internal states. This could be controlled by a user defined switch.
      int numStates = geomState_->getNumElemStates(domain_->getElementSet()[iElem]->getGlNum());
      if(numStates > 0) {
        double *states = geomState_->getElemState(domain_->getElementSet()[iElem]->getGlNum());
        for (int j = 0; j < numStates; ++j) states[j] = 0;
        domain_->getElementSet()[iElem]->initStates(states);
      }
    }

    delete [] nodes;
  }
  filePrint(stderr,"\r %4.2f%% complete\n", 100.);
}

template<typename MatrixBufferType, typename SizeType>
void
ElementSamplingDriver<MatrixBufferType,SizeType>::solve() {

  AllOps<double> allOps;
  preProcessGlobal(allOps);

  std::vector<int> sampleElemIds, packedToInput;
  std::vector<std::map<int, double> > weights(domain_->solInfo().readInROBorModes.size());
  makePackedToInput(packedToInput);

  for(int j=0; j<domain_->solInfo().readInROBorModes.size(); ++j) {

    preProcessLocal(allOps, j);

    // Training target (solver_.rhsBuffer) is the sum of elementary contributions
#ifdef PRINT_ESTIMERS
    double t2 = getTime();
#endif
    for(int i=0; i<solver_.equationCount(); ++i) solver_.rhsBuffer()[i] = 0.0;
    assembleTrainingData(podBasis_, podBasis_.vectorCount(), displac_, veloc_, accel_, j);
#ifdef PRINT_ESTIMERS
    fprintf(stderr, "time for assembleTrainingData = %f\n", (getTime()-t2)/1000.0);
#endif

    Vector solution;
    computeSolution(solution, domain_->solInfo().tolPodRom);

    postProcessLocal(solution, packedToInput, j, sampleElemIds, weights[j]);
  }

  postProcessGlobal(sampleElemIds, weights);
}

template<typename MatrixBufferType, typename SizeType>
void
ElementSamplingDriver<MatrixBufferType,SizeType>::computeSolution(Vector &solution, double relativeTolerance, bool verboseFlag) {

  solver_.relativeToleranceIs(relativeTolerance);
  solver_.verboseFlagIs(verboseFlag);
  solver_.scalingFlagIs(domain->solInfo().useScalingSpnnls);
  solver_.centerFlagIs(domain->solInfo().useCenterSpnnls);
  solver_.reverseFlagIs(domain->solInfo().useReverseOrder);
  solver_.projectFlagIs(domain->solInfo().projectSolution);
  solver_.positivityIs(domain->solInfo().positiveElements);
  solver_.solverTypeIs(domain->solInfo().solverTypeSpnnls);
  solver_.maxSizeRatioIs(domain->solInfo().maxSizeSpnnls);
  solver_.maxIterRatioIs(domain->solInfo().maxIterSpnnls);
  solver_.maxNumElemsIs(domain->solInfo().maxElemSpnnls);
  try {
    solver_.solve();
  }
  catch(std::runtime_error& e) {
    std::cerr << " *** WARNING: " << e.what() << std::endl;
  }

  if(verboseFlag) {
    std::cout << "Primal solution:";
    for (int elemRank = 0; elemRank != elementCount(); ++elemRank) {
      std::cout << " " << solver_.solutionEntry(elemRank);
    }
    std::cout << "\n";

    StackVector trainingTarget(solver_.rhsBuffer(), solver_.equationCount());
    std::cout << "Error magnitude / Absolute tolerance = " << solver_.errorMagnitude() << " / " << solver_.relativeTolerance() * trainingTarget.norm() << "\n";
    std::cout << "1-norm of primal solution = " << std::accumulate(solver_.solutionBuffer(), solver_.solutionBuffer() + solver_.unknownCount(), 0.0) << "\n";
  }

  // Read solution
  solution.initialize(elementCount());
  std::copy(solver_.solutionBuffer(), solver_.solutionBuffer() + elementCount(), solution.data());
}

template<typename MatrixBufferType, typename SizeType>
void
ElementSamplingDriver<MatrixBufferType,SizeType>::postProcess(Vector &solution, bool verboseFlag) {

  std::vector<int> sampleElemIds, packedToInput;
  std::vector<std::map<int, double> > weights(1);

  makePackedToInput(packedToInput);
  postProcessLocal(solution, packedToInput, 0, sampleElemIds, weights[0], verboseFlag);
  postProcessGlobal(sampleElemIds, weights, verboseFlag);
}

template<typename MatrixBufferType, typename SizeType>
void
ElementSamplingDriver<MatrixBufferType,SizeType>::makePackedToInput(std::vector<int> &packedToInput)
{
  packedToInput.resize(elementCount());
  Elemset &inputElemSet = *(geoSource->getElemSet());
  for (int iElem = 0, iElemEnd = inputElemSet.size(); iElem != iElemEnd; ++iElem) {
    Element *elem = inputElemSet[iElem];
    if (elem) {
      const int iPackElem = domain_->glToPackElem(iElem);
      if(iPackElem >= 0) {
        assert(iPackElem < packedToInput.size());
        packedToInput[iPackElem] = iElem;
      }
    }
  }
}

template<typename MatrixBufferType, typename SizeType>
void
ElementSamplingDriver<MatrixBufferType,SizeType>::postProcessLocal(Vector &solution, std::vector<int> packedToInput, int j,
                                                                   std::vector<int> &sampleElemIds, std::map<int, double> &weights,
                                                                   bool verboseFlag)
{
  std::set<int> sampleElemRanks;
  {
    for (int iElem = 0; iElem != elementCount(); ++iElem) {
      if (solution[iElem] > 0.0) {
        sampleElemRanks.insert(sampleElemRanks.end(), iElem);
      }
    }
  }

  sampleElemIds.reserve(sampleElemRanks.size());
  for (std::set<int>::const_iterator it = sampleElemRanks.begin(), it_end = sampleElemRanks.end(); it != it_end; ++it) {
    const int elemRank = packedToInput[*it];
    weights.insert(std::make_pair(elemRank, solution[*it]));
    if(domain_->solInfo().readInROBorModes.size() == 1 || std::find(sampleElemIds.begin(), sampleElemIds.end(), elemRank) == sampleElemIds.end())
      sampleElemIds.push_back(elemRank);
  }

  outputFullWeights(solution, packedToInput, ((domain_->solInfo().readInROBorModes.size() == 1) ? -1 : j));
}

template<typename MatrixBufferType, typename SizeType>
void
ElementSamplingDriver<MatrixBufferType,SizeType>::postProcessGlobal(std::vector<int> &sampleElemIds, std::vector<std::map<int, double> > &weights,
                                                                    bool verboseFlag)
{
  const FileNameInfo fileInfo;

  if(domain_->solInfo().localBasisSize.size() > 1) {
    int globalBasisSize = std::accumulate(domain_->solInfo().localBasisSize.begin(), domain_->solInfo().localBasisSize.end(), 0);
    podBasis_.localBasisIs(0, globalBasisSize);
  }

  // compute the reduced forces
  Vector forceFull(SingleDomainDynamic::solVecInfo(),0.0);
  // 1) gravity
  Vector gravForceRed(podBasis_.vectorCount());
  if(domain->gravityFlag()) {
    domain->addGravityForce(forceFull);
    reduce(podBasis_, forceFull, gravForceRed);
  }
  // 2) constant force or constant part of time-dependent forces (default loadset only) TODO add support for multiple loadsets
  domain->computeUnamplifiedExtForce(forceFull, 0);
  Vector constForceRed(podBasis_.vectorCount());
  bool reduce_f = (forceFull.norm() != 0);
  if(reduce_f) reduce(podBasis_, forceFull, constForceRed);

  // 3) set LMPC
  int numLMPC = domain->getNumLMPC(); 

  // compute the reduced initial conditions
  Vector d0Full(SingleDomainDynamic::solVecInfo()),
         v0Full(SingleDomainDynamic::solVecInfo());
  Vector tmp(SingleDomainDynamic::solVecInfo());
  SysState<Vector> inState(d0Full, v0Full, tmp, tmp);
  SingleDomainDynamic::getInitState(inState);
  Vector d0Red(podBasis_.vectorCount()),
         v0Red(podBasis_.vectorCount());
  bool reduce_idis = (d0Full.norm() != 0),
       reduce_ivel = (v0Full.norm() != 0);
#ifdef USE_EIGEN3
  GenEiSparseMatrix<double,Eigen::SimplicialLLT<Eigen::SparseMatrix<double>,Eigen::Upper> > *M;
#endif
  if(domain->solInfo().useMassNormalizedBasis || domain->solInfo().newmarkBeta == 0 || domain_->solInfo().readInROBorModes.size() > 1) {
    AllOps<double> allOps;
    if(reduce_idis || reduce_ivel || domain_->solInfo().readInROBorModes.size() > 1) { 
#ifdef USE_EIGEN3
      allOps.M = M = domain_->constructEiSparseMatrix<double,Eigen::SimplicialLLT<Eigen::SparseMatrix<double>,Eigen::Upper> >();
#else
      allOps.M = domain->constructDBSparseMatrix<double>();
#endif
      domain->makeSparseOps(allOps, 0, 0, 0);
    }
    if(reduce_idis) {
      allOps.M->mult(d0Full, tmp);
      reduce(podBasis_, tmp, d0Red);
    }
    if(reduce_ivel) {
      allOps.M->mult(v0Full, tmp);
      reduce(podBasis_, tmp, v0Red);
    }
  }
  else {
    if(reduce_idis) reduce(podBasis_, d0Full, d0Red);
    if(reduce_ivel) reduce(podBasis_, v0Full, v0Red);
  }

  // output the reduced mesh
  Elemset &inputElemSet = *(geoSource->getElemSet());
  std::auto_ptr<Connectivity> elemToNode(new Connectivity(&inputElemSet));
  const MeshRenumbering meshRenumbering(sampleElemIds.begin(), sampleElemIds.end(), *elemToNode, verboseFlag);
  const MeshDesc reducedMesh(domain_, geoSource, meshRenumbering, weights);
  if(domain_->solInfo().localBasisSize.size() == 1)
    outputMeshFile(fileInfo, reducedMesh, podBasis_.vectorCount());
  else
    outputMeshFile(fileInfo, reducedMesh, domain_->solInfo().localBasisSize);

  // output the reduced forces
  std::ofstream meshOut(getMeshFilename(fileInfo).c_str(), std::ios_base::app);
  if(domain->solInfo().reduceFollower) meshOut << "*\nEXTFOL\n";
  if(domain->gravityFlag()) {
    meshOut << "*\nFORCES -1\nMODAL\n"; // note: gravity forces are put in loadset -1 so that MFTT (if present) will not be applied
    meshOut.precision(std::numeric_limits<double>::digits10+1);
    for(int i=0; i<podBasis_.vectorCount(); ++i)
      meshOut << i+1 << " " << gravForceRed[i] << std::endl;
  }
  if(reduce_f) {
    meshOut << "*\nFORCES\nMODAL\n";
    meshOut.precision(std::numeric_limits<double>::digits10+1);
    for(int i=0; i<podBasis_.vectorCount(); ++i) 
      meshOut << i+1 << " " << constForceRed[i] << std::endl;
  }

#ifdef USE_EIGEN3
  if(numLMPC > 0) {

    if(strcmp(domain_->solInfo().readInDualROB,"") != 0) {
      VecBasis dualProjectionBasis_;
      VecNodeDof1Conversion vecNodeDof1Conversion(numLMPC);
      // Load dual projection basis    
      std::string fileName = BasisFileId(fileInfo, BasisId::DUALSTATE, BasisId::POD);
      BasisInputStream<1> dualProjectionBasisInput(fileName, vecNodeDof1Conversion);
      const int dualProjectionSubspaceSize = domain_->solInfo().maxSizeDualBasis ?
                                             std::min(domain_->solInfo().maxSizeDualBasis, dualProjectionBasisInput.size()) :
                                             dualProjectionBasisInput.size();
  
      readVectors(dualProjectionBasisInput, dualProjectionBasis_, dualProjectionSubspaceSize);
  
      filePrint(stderr, " ... Dual Proj. Subspace Dim. = %-3d ...\n", dualProjectionBasis_.vectorCount());
      meshOut << "*\nLMPC\nMODAL " << dualProjectionBasis_.vectorCount() << std::endl;

      std::vector<Eigen::Triplet<double> > tripletList;
      Eigen::SparseMatrix<double> C(numLMPC, domain_->getCDSA()->size());
      Eigen::Matrix<double,Eigen::Dynamic,1> g(numLMPC);
      LMPCons** lmpc(domain_->getLMPC()->data());
      
      //construct constraint matrix and right hand side
      for(int i=0; i<numLMPC; ++i) {
        for(int j=0; j<lmpc[i]->nterms; ++j) {
          int cdof = domain_->getCDSA()->locate(lmpc[i]->terms[j].nnum, 1 << lmpc[i]->terms[j].dofnum);
          if(cdof > -1) {
            tripletList.push_back(Eigen::Triplet<double>(i, cdof, double(lmpc[i]->terms[j].coef.r_value)));
          }
        }
        g[i] = lmpc[i]->rhs.r_value;
      } 
 
      C.setFromTriplets(tripletList.begin(), tripletList.end());
      const Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,Eigen::ColMajor> &V = podBasis_.basis(), &W = dualProjectionBasis_.basis();
      Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> reducedConstraintMatrix_ = W.transpose()*C*V;
      Eigen::Matrix<double,Eigen::Dynamic,1> reducedConstraintRhs0_ = W.transpose()*g;
      // output reduced Constraints to reduced mesh file
      for(int col=0; col<podBasis_.vectorCount(); ++col) {
         meshOut << reducedConstraintMatrix_.col(col) << std::endl;
      }
      meshOut << reducedConstraintRhs0_ << std::endl;

    } else {
      filePrint(stderr, "... No Filename given for Precomputation of reduced Constraints ...\n");
    }

  } 
#endif

  // output the reduced initial conditions
  if(reduce_idis) {
    meshOut << "*\nIDISPLACEMENTS\nMODAL\n";
    meshOut.precision(std::numeric_limits<double>::digits10+1);
    for(int i=0; i<podBasis_.vectorCount(); ++i)
      meshOut << i+1 << " " << d0Red[i] << std::endl;
  }
  if(reduce_ivel) {
    meshOut << "*\nIVELOCITIES\nMODAL\n";
    meshOut.precision(std::numeric_limits<double>::digits10+1);
    for(int i=0; i<podBasis_.vectorCount(); ++i)
      meshOut << i+1 << " " << v0Red[i] << std::endl;
  }

#ifdef USE_EIGEN3
  // pre-computation required for local bases method
  if(domain_->solInfo().readInROBorModes.size() > 1 && 
     (domain->solInfo().readInLocalBasesCent.size() == domain_->solInfo().readInROBorModes.size())) {
    // read cluster centroids
    const int n = podBasis_.vectorSize();
    Eigen::MatrixXd uc(n, domain->solInfo().readInLocalBasesCent.size());
    const VecNodeDof6Conversion vecNodeDof6Conversion(*domain_->getCDSA());
    if(domain->solInfo().readInLocalBasesCent.size() == domain_->solInfo().readInROBorModes.size()) {
      for(int i=0; i<domain->solInfo().readInLocalBasesCent.size(); ++i) {
        BasisInputStream<6> in(domain->solInfo().readInLocalBasesCent[i], vecNodeDof6Conversion);
        in >> uc.col(i).data();
      }
    }
 
    // compute and output auxiliary quantities
    meshOut << "*\nLOCROB\n";
    for(int i=0; i<domain_->solInfo().readInROBorModes.size(); ++i) {
      int startColi = std::accumulate(domain_->solInfo().localBasisSize.begin(), domain_->solInfo().localBasisSize.begin()+i, 0);
      int blockColsi = domain_->solInfo().localBasisSize[i];
      for(int j=i+1; j<domain_->solInfo().readInROBorModes.size(); ++j) {
        int startColj = std::accumulate(domain_->solInfo().localBasisSize.begin(), domain_->solInfo().localBasisSize.begin()+j, 0);
        int blockColsj = domain_->solInfo().localBasisSize[j];
        Eigen::MatrixXd VtVij;
        std::string fileName = domain->solInfo().reducedMeshFile;
        std::ostringstream ss;
        ss << ".auxiliary.cluster" << i+1 << "." << j+1;
        fileName.append(ss.str());
        meshOut << "auxi " << i+1 << " " << j+1 << " \"" << fileName << "\"\n";
        filePrint(stderr," ... Writing local bases auxiliary quantities to file %s ...\n", fileName.c_str());
        std::ofstream matrixOut(fileName);
        // matrix: Vi.transpose()*Vj (or Vi.transpose()*M*Vj)
        if(domain_->solInfo().newmarkBeta == 0 || domain_->solInfo().useMassNormalizedBasis) {
          VtVij = podBasis_.basis().block(0,startColi,n,blockColsi).transpose()*
                  (M->getEigenSparse().selfadjointView<Eigen::Upper>()*podBasis_.basis().block(0,startColj,n,blockColsj)); 
        }
        else {
           VtVij = podBasis_.basis().block(0,startColi,n,blockColsi).transpose()*
                   podBasis_.basis().block(0,startColj,n,blockColsj);
        }
        matrixOut << std::setprecision(16) << VtVij << std::endl;
        // scalar: ||uci|| - ||ucj||
        double d = uc.col(i).squaredNorm() - uc.col(j).squaredNorm();
        matrixOut << d << std::endl;
        // vector: 2*V.transpose()*(ucj-uci)
        Eigen::VectorXd w = 2*podBasis_.basis().transpose()*(uc.col(j)-uc.col(i));
        matrixOut << w.transpose() << std::endl;

        matrixOut.close();
      }
    }
  }

  // build and output compressed basis
  DofSetArray reduced_dsa(reducedMesh.nodes().size(), const_cast<Elemset&>(reducedMesh.elements()));
  int num_bc = reducedMesh.dirichletBConds().size();
  BCond *bc = (num_bc > 0) ? const_cast<BCond*>(&reducedMesh.dirichletBConds()[0]) : NULL;
  ConstrainedDSA reduced_cdsa(reduced_dsa, num_bc, bc);
  podBasis_.makeSparseBasis(meshRenumbering.reducedNodeIds(), domain_->getCDSA(), &reduced_cdsa);
  {
    VecNodeDof6Conversion converter(reduced_cdsa);
    for(int j=0; j<domain_->solInfo().readInROBorModes.size(); ++j) {
      std::string filename = getMeshFilename(fileInfo).c_str();
      filename.append(".compressed.basis");
      if(domain_->solInfo().readInROBorModes.size() != 1) {
        std::ostringstream ss;
        ss << j+1;
        filename.append(ss.str());
      }
      if(domain_->solInfo().newmarkBeta == 0 || domain_->solInfo().useMassNormalizedBasis) filename.append(".normalized");
      filePrint(stderr," ... Writing compressed basis to file %s ...\n", filename.c_str());
      BasisOutputStream<6> output(filename, converter, false);

      int startCol = (domain_->solInfo().readInROBorModes.size() == 1) ? 0 :
                      std::accumulate(domain_->solInfo().localBasisSize.begin(), domain_->solInfo().localBasisSize.begin()+j, 0);
      int blockCols = (domain_->solInfo().readInROBorModes.size() == 1) ? podBasis_.vectorCount() : domain_->solInfo().localBasisSize[j];
      for (int iVec = startCol; iVec < startCol+blockCols; ++iVec) {
        output << podBasis_.compressedBasis().col(iVec);
      }
    }
  }
#endif
}

template<typename MatrixBufferType, typename SizeType>
void
ElementSamplingDriver<MatrixBufferType,SizeType>::preProcess()
{
  AllOps<double> allOps;
  preProcessGlobal(allOps);
  preProcessLocal(allOps, 0);
}

template<typename MatrixBufferType, typename SizeType>
void
ElementSamplingDriver<MatrixBufferType,SizeType>::preProcessGlobal(AllOps<double>& allOps)
{
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
  if(domain_->solInfo().newmarkBeta == 0) {
    domain_->assembleNodalInertiaTensors(melArray_);
  }

  // Assemble mass matrix if necessary
  if(domain_->solInfo().useMassOrthogonalProjection) {
    if(geoSource->getMRatio() != 0) {
      allOps.M = domain_->constructDBSparseMatrix<double>();
    }
    else {
      allOps.M = new DiagMatrix(domain->getCDSA());
    }
    domain_->makeSparseOps<double>(allOps, 0.0, 1.0, 0.0);
  }
}

template<typename MatrixBufferType, typename SizeType>
void
ElementSamplingDriver<MatrixBufferType,SizeType>::preProcessLocal(AllOps<double> &allOps, int j)
{
  const FileNameInfo fileInfo;
  
  // Read order reduction data
  const VecNodeDof6Conversion vecDofConversion(*domain_->getCDSA());
  assert(vectorSize() == vecDofConversion.vectorSize());

  // Read in basis to be used for the projection:
  // (a) if a mass-orthogonal projection is to be done, then the mass-normalized basis will be read
  // (b) if and orthogonal projection is to be done, then the identity-normalized basis will be read
  {
    std::string fileName = BasisFileId(fileInfo, BasisId::STATE, BasisId::POD, j);
    if(domain_->solInfo().useMassOrthogonalProjection) fileName.append(".normalized");
    BasisInputStream<6> in(fileName, vecDofConversion);
    if(domain_->solInfo().readInROBorModes.size() == 1) {
      const int podSizeMax = domain_->solInfo().maxSizePodRom;
      if(podSizeMax != 0) {
        readVectors(in, podBasis_, podSizeMax);
      } else {
        readVectors(in, podBasis_);
      }
    }
    else {
      int globalBasisSize = std::accumulate(domain_->solInfo().localBasisSize.begin(), domain_->solInfo().localBasisSize.end(), 0);
      int startCol = std::accumulate(domain_->solInfo().localBasisSize.begin(), domain_->solInfo().localBasisSize.begin()+j, 0);
      int blockCols = domain_->solInfo().localBasisSize[j];
      readVectors(in, podBasis_, globalBasisSize, blockCols, startCol);
      podBasis_.localBasisIs(startCol, blockCols);
    }
  }

  const int podVectorCount = podBasis_.vectorCount();

  // Read some displacement snapshots from one or more files and project them on to the basis
  readAndProjectSnapshots(BasisId::STATE, vectorSize(), podBasis_, vecDofConversion,
                          snapshotCounts_, timeStamps_, displac_, allOps.M,
                          ((domain_->solInfo().readInROBorModes.size() == 1) ? -1 : j));

  // Optionally, read some velocity snapshots and project them on to the reduced order basis
  if(!domain_->solInfo().velocPodRomFile.empty()) {
    std::vector<double> velTimeStamps;
    std::vector<int> velSnapshotCounts;
    veloc_ = new VecBasis;
    readAndProjectSnapshots(BasisId::VELOCITY, vectorSize(), podBasis_, vecDofConversion,
                            velSnapshotCounts, velTimeStamps, *veloc_, allOps.M,
                            ((domain_->solInfo().readInROBorModes.size() == 1) ? -1 : j));
    if(velSnapshotCounts != snapshotCounts_) std::cerr << " *** WARNING: inconsistent velocity snapshots\n";
  }

  // Optionally, read some acceleration snapshots and project them on to the reduced order basis
  if(!domain_->solInfo().accelPodRomFile.empty()) {
    std::vector<double> accTimeStamps;
    std::vector<int> accSnapshotCounts;
    accel_ = new VecBasis;
    readAndProjectSnapshots(BasisId::ACCELERATION, vectorSize(), podBasis_, vecDofConversion,
                            accSnapshotCounts, accTimeStamps, *accel_, allOps.M,
                            ((domain_->solInfo().readInROBorModes.size() == 1) ? -1 : j));
    if(accSnapshotCounts != snapshotCounts_) std::cerr << " *** WARNING: inconsistent acceleration snapshots\n";
  }
  
  // Read in mass-normalized basis if necessary
  if((domain_->solInfo().newmarkBeta == 0 || domain_->solInfo().useMassNormalizedBasis)
     && !domain_->solInfo().useMassOrthogonalProjection) {
    std::string fileName = BasisFileId(fileInfo, BasisId::STATE, BasisId::POD,j);
    fileName.append(".normalized");
    BasisInputStream<6> in(fileName, vecDofConversion);
    if(domain_->solInfo().readInROBorModes.size() == 1) {
      const int podSizeMax = domain_->solInfo().maxSizePodRom;
      if(podSizeMax != 0) {
        readVectors(in, podBasis_, podSizeMax);
      } else {
        readVectors(in, podBasis_);
      }
    }
    else {
      int globalBasisSize = std::accumulate(domain_->solInfo().localBasisSize.begin(), domain_->solInfo().localBasisSize.end(), 0);
      int startCol = std::accumulate(domain_->solInfo().localBasisSize.begin(), domain_->solInfo().localBasisSize.begin()+j, 0);
      int blockCols = domain_->solInfo().localBasisSize[j];
      readVectors(in, podBasis_, globalBasisSize, blockCols, startCol);
    }
  }

  const int snapshotCount = (j==-1) ? std::accumulate(snapshotCounts_.begin(), snapshotCounts_.end(), 0) : snapshotCounts_[j];
  solver_.problemSizeIs(podVectorCount*snapshotCount, elementCount());
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
  return new Rom::ElementSamplingDriver<std::vector<double>,size_t>(d);
}
