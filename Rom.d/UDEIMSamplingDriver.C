#include "UDEIMSamplingDriver.h"

#include "ElementSamplingDriver.h"

#include "SvdOrthogonalization.h"
#include "VecNodeDof6Conversion.h"
#include "NodalRestrictionMapping.h"
#include "ConnectivityUtils.h"
#include "VecNodeDof6Map.h"
#include "BasisFileStream.h"
#include "FileNameInfo.h"
#include "SimpleBuffer.h"

#include "RenumberingUtils.h"
#include "MeshDesc.h"

#include "DistrBasisFile.h"
#include <Timers.d/StaticTimers.h>

#include <Driver.d/Domain.h>
#include <Driver.d/SysState.h>
#include <Driver.d/GeoSource.h>
#include <Utils.d/dofset.h>
#include <Utils.d/DistHelper.h>

#ifdef USE_EIGEN3
#include <Eigen/Dense>
#endif
#include <cmath>
#include <utility>
#include <algorithm>
#include <numeric>

namespace Rom {

UDEIMSamplingDriver::UDEIMSamplingDriver(Domain *domain) :
  SingleDomainDynamic(domain),
  converter(NULL)
{}

void
UDEIMSamplingDriver::solve() {

  SingleDomainDynamic::preProcess();
  converter = new VecNodeDof6Conversion(*domain->getCDSA());

  const int podSizeMax = domain->solInfo().maxSizePodRom; 
  bool normalized = true;

  if(domain->solInfo().computeForceSnap){
    readInBasis(podBasis_, BasisId::STATE, BasisId::POD, podSizeMax);
    writeProjForceSnap();
  } else if(domain->solInfo().orthogForceSnap) {
    int forcePodSizeMax = domain->solInfo().forcePodSize;
    VecBasis forceBasis;
    readInBasis(forceBasis, BasisId::FORCE, BasisId::SNAPSHOTS,forcePodSizeMax);
    std::vector<double> SVs;
    OrthoForceSnap(forceBasis,SVs);
    writeBasisToFile(forceBasis, SVs, BasisId::FORCE, BasisId::ROB);
  } else {
    readInBasis(podBasis_, BasisId::STATE, BasisId::POD, podSizeMax, normalized);
    std::vector<int> maskIndices;
    VecBasis forceBasis;
    computeInterpIndices(forceBasis, maskIndices);
    computeAndWriteUDEIMBasis(forceBasis, maskIndices);
    writeSampledMesh(maskIndices);
    }

}

void
UDEIMSamplingDriver::readInBasis(VecBasis &podBasis, BasisId::Type type, BasisId::Level level, int podSizeMax, bool normalized)
{
 FileNameInfo fileInfo;
 std::string fileName = BasisFileId(fileInfo, type, level);
 if(normalized) fileName.append(".normalized");
 BasisInputStream in( fileName, *converter);
 if (podSizeMax != 0) {
   std::cout << "reading in " << podSizeMax << " vectors from " << fileName.c_str() << std::endl;
   readVectors(in, podBasis, podSizeMax);
 } else {
   std::cout << "reading in all vectors from " << fileName.c_str() << std::endl;
   readVectors(in, podBasis);
 }
}

template <typename Scalar>
void
UDEIMSamplingDriver::writeBasisToFile(const VecBasis &OutputBasis, std::vector<Scalar> singularValue, BasisId::Type type, BasisId::Level level)
{
 FileNameInfo fileInfo;
 std::string fileName = BasisFileId(fileInfo, type, level);
 BasisOutputStream outputNormalized(fileName, *converter, false);
 filePrint(stderr, " ... Writing basis to file %s ...\n", fileName.c_str());
 for (int iVec = 0; iVec < OutputBasis.vectorCount(); ++iVec) {
   if(singularValue.size() > 0)
     outputNormalized << std::make_pair(double(singularValue[iVec]), OutputBasis[iVec]);
   else
     outputNormalized << OutputBasis[iVec];
 }
}

void
UDEIMSamplingDriver::writeProjForceSnap() 
{
  VecBasis forceBasis;
 
  //First read in state snapshots 
  VecBasis displac;
  std::vector<double> timeStamps;
  std::vector<int> snapshotCounts;
  readAndProjectSnapshots(BasisId::STATE, converter->vectorSize(), podBasis_, converter,
                          snapshotCounts, timeStamps, displac);

  VecBasis *veloc_;
  if(!domain->solInfo().velocPodRomFile.empty()) {
   std::vector<double> velTimeStamps;
   std::vector<int> velSnapshotCounts;
   veloc_ = new VecBasis;
   readAndProjectSnapshots(BasisId::VELOCITY, converter->vectorSize(), podBasis_, converter,
                           velSnapshotCounts, velTimeStamps, *veloc_);
   if(velSnapshotCounts != snapshotCounts) std::cerr << " *** WARNING: inconsistent velocity snapshots\n";
  } else { veloc_ = NULL;}

  VecBasis *accel_;
  if(!domain->solInfo().accelPodRomFile.empty()) {
   std::vector<double> accTimeStamps;
   std::vector<int> accSnapshotCounts;
   accel_ = new VecBasis;
   readAndProjectSnapshots(BasisId::ACCELERATION, converter->vectorSize(), podBasis_, converter,
                           accSnapshotCounts, accTimeStamps, *accel_);
   if(accSnapshotCounts != snapshotCounts) std::cerr << " *** WARNING: inconsistent acceleration snapshots\n";
  } else { accel_ = NULL; } 

  //Now build forcevectors
  buildForceArray(forceBasis, displac, veloc_, accel_, timeStamps, snapshotCounts);
  writeBasisToFile(forceBasis, timeStamps, BasisId::FORCE, BasisId::SNAPSHOTS);
}

void
UDEIMSamplingDriver::computeInterpIndices(VecBasis &forceBasis, std::vector<int> &maskIndices) {
  //member function for determining the sampling indicies for UDEIM
  //result is (U*(P^T*U)^-1) where U is the left singular vectors of force snapshots and P
  //is the column selection matrix
  
  int forcePodSizeMax = domain->solInfo().forcePodSize;
  readInBasis(forceBasis, BasisId::FORCE, BasisId::SNAPSHOTS,forcePodSizeMax);  
#ifdef USE_EIGEN3
  Eigen::Map< Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> > forceMatrix(forceBasis.data(),forceBasis.vectorSize(),forceBasis.vectorCount());

  int maxCoeffSlot;

  {
   Eigen::Matrix<double,Eigen::Dynamic,1> firstCol(forceMatrix.rows());
   firstCol = forceMatrix.col(0);
   //take absolute value of first column
   for(int row = 0; row != firstCol.rows(); row++)
     if(firstCol(row) < 0.0)
       firstCol(row) = -1.0*firstCol(row);
   //find maximum component
   firstCol.maxCoeff(&maxCoeffSlot);
   maskIndices.push_back(maxCoeffSlot);
  }

  //start loop to compute mask indicies 
  for(int i = 1; i < forceMatrix.cols(); ++i){ //loop starts at 1 i.e. the 2nd column
    filePrint(stderr,"\r %4.2f%% complete", double(i)/double(forceBasis.vectorCount())*100.);

    //allocate space for P^T*U and P^T*u_i
    Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> Umasked(maskIndices.size(),maskIndices.size());
    Eigen::Matrix<double,Eigen::Dynamic,1>              u_i_masked(maskIndices.size());

    for(int j = 0; j < maskIndices.size(); ++j) {//select proper rows and columns of force basis
      Umasked.row(j) = forceMatrix.block(maskIndices[j],0,1,maskIndices.size()); //(P^T*U) is square
      u_i_masked(j) = forceMatrix(maskIndices[j],i);//mask next column over
    }

    Eigen::FullPivLU< Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> > luOfUmasked(Umasked); //invert row reduced basis

    Eigen::Matrix<double,Eigen::Dynamic,1> residual(forceBasis.vectorSize());

    if(luOfUmasked.isInvertible())
      residual = forceMatrix.col(i) - forceMatrix.leftCols(maskIndices.size())*luOfUmasked.inverse()*u_i_masked;
    else
      throw std::runtime_error("... Matrix Not Invertible ...");

    //take absolute value of residual vector
    for(int row = 0; row != residual.rows(); row++)
      if(residual(row) < 0.0)
        residual(row) = -1.0*residual(row);

    residual.maxCoeff(&maxCoeffSlot); //find indice of maximum component
    maskIndices.push_back(maxCoeffSlot);
  }
  filePrint(stderr,"\r %4.2f%% complete\n", 100.);

  Eigen::Map< Eigen::Matrix<int,Eigen::Dynamic,Eigen::Dynamic> > indSol(maskIndices.data(),maskIndices.size(),1);
  std::cout << "selected indices:" << indSol.transpose() << std::endl;
#endif
}

void
UDEIMSamplingDriver::computeAndWriteUDEIMBasis(VecBasis &forceBasis, std::vector<int> &maskIndices)
{
  //member function for computing and writing the UDEIM basis P*(P^T*U)^-T*U^T*V
  //where V is the mass-orthogonal POD basis, U is the left Singular vectors of the force snapshots
  //and P column selection matrix derived form the sampled indicies computed above
  int maxDeimBasisSize = domain->solInfo().maxDeimBasisSize;
  if(maxDeimBasisSize == 0)
    maxDeimBasisSize = maskIndices.size();

  VecBasis deimBasis(podBasis_.vectorCount(),podBasis_.vectorInfo());
#ifdef USE_EIGEN3
  Eigen::Map<Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> > podMap(podBasis_.data(),podBasis_.size(), podBasis_.vectorCount());
  Eigen::Map<Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> > forceMap(forceBasis.data(),forceBasis.size(), forceBasis.vectorCount());
  Eigen::Map<Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> > deimMap(deimBasis.data(),deimBasis.size(), deimBasis.vectorCount());

  Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> compressedDBTranspose(podBasis_.numVectors(),maskIndices.size());
  Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> rowReducedFM(maskIndices.size(),maxDeimBasisSize);

  //initialize (P^T*U)
  for(int row = 0; row != maskIndices.size(); row++)
    rowReducedFM.row(row) = forceMap.block(maskIndices[row],0,1,maxDeimBasisSize);

  //compute W^T = V^T*U*(P^T*U)^-1
  Eigen::JacobiSVD< Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> > SVDOfUmasked(rowReducedFM,Eigen::ComputeThinU | Eigen::ComputeThinV);
  Eigen::Matrix<double,Eigen::Dynamic,1> invSVs(SVDOfUmasked.nonzeroSingularValues());

  for(int i = 0; i != invSVs.rows(); i++) invSVs(i) = 1.0/SVDOfUmasked.singularValues()(i);

  compressedDBTranspose = podMap.transpose()*forceMap.leftCols(maxDeimBasisSize)*SVDOfUmasked.matrixV()*invSVs.asDiagonal()*SVDOfUmasked.matrixU().transpose();
  //we are computing the transpose of the basis

  std::cout << "compressed Basis" << std::endl;
  std::cout << compressedDBTranspose.transpose() << std::endl;

  //initialize deim basis to all zeros
  for(int i = 0; i != deimMap.rows(); i++)
    for(int j = 0; j != deimMap.cols(); j++)
      deimMap(i,j) = 0.0;

  //expand rows W^T*P^T
  std::cout << " num of rows = " << compressedDBTranspose.rows() << " num of cols = " << compressedDBTranspose.cols() << std::endl;
  for(int row = 0; row != maskIndices.size(); row++) {
    deimMap.row(maskIndices[row]) = compressedDBTranspose.col(row);//use .col member to get rows of transposed basis
  }

/*  std::cout << "deimMap" << std::endl;
  std::cout << deimMap << std::endl;

  std::cout << "podMap" << std::endl;
  std::cout << podMap << std::endl;*/

  std::vector<double> dummySVs; 
  writeBasisToFile(deimBasis, dummySVs, BasisId::FORCE, BasisId::ROB);
#endif
}

void
UDEIMSamplingDriver::writeSampledMesh(std::vector<int> &maskIndices) {

  VecNodeDof6Map nodeDofMap(*domain->getCDSA());

  //get nodes & dofs belonging to sampled indices (global numbering)
  std::set<int> selectedNodeSet; //To ensure unicity, since several locations (i.e. vector indices) can correspond to one node
  std::vector<std::pair<int,int> > compressedNodeKey; //need separate vector container for repeated nodes with different selected dofs
  for (std::vector<int>::const_iterator it = maskIndices.begin(); it != maskIndices.end(); ++it) {
    selectedNodeSet.insert(nodeDofMap.nodeDof(*it).nodeRank);
    int selectedNode = nodeDofMap.nodeDof(*it).nodeRank;
    int selectedNodeDof = std::log10(nodeDofMap.nodeDof(*it).dofId)/std::log10(2.0); //this function call is returning 2^dof for some reason, call log to normalize
    compressedNodeKey.push_back(std::make_pair(selectedNode,selectedNodeDof));
  }

//  std::sort(compressedNodeKey.begin(),compressedNodeKey.end());

  //print nodes to screen
  filePrint(stderr,"selected Nodes:");
  for(std::set<int>::iterator it = selectedNodeSet.begin(); it != selectedNodeSet.end(); it++)
    filePrint(stderr," %d", *it);

  filePrint(stderr,"\n");

  filePrint(stderr,"number of selectedNodeSet = %d, number of Indices = %d, number of dofs = %d\n", selectedNodeSet.size(),maskIndices.size(), compressedNodeKey.size());

  {
   DofSetArray *cdsa = domain->getCDSA();
   for(std::vector<std::pair<int,int> >::iterator it = compressedNodeKey.begin(); it != compressedNodeKey.end(); it++)
     std::cout << "node: " << it->first+1 << " dof1 = " << cdsa->firstdof(it->first)+1 << " slot = " << it->second << std::endl;
  }

  // compute the reduced forces (constant only)
  Vector constForceFull(solVecInfo());
  getConstForce(constForceFull);
  Vector constForceRed(podBasis_.vectorCount());
  reduce(podBasis_, constForceFull,  constForceRed);

  // Determine mapping between elements and nodes
  std::auto_ptr<Connectivity> elemToNode(new Connectivity(geoSource->getElemSet()));
  std::auto_ptr<Connectivity> nodeToElem(elemToNode->reverse());

  //fill element container with all elements 
  std::vector<int> packedToInput(elementCount());
  {
   Elemset &inputElemSet = *(geoSource->getElemSet());
   for (int iElem = 0, iElemEnd = inputElemSet.size(); iElem != iElemEnd; ++iElem) {
     Element *elem = inputElemSet[iElem];
     if (elem) {
       const int iPackElem = domain->glToPackElem(iElem);
       if(iPackElem >= 0) {
         assert(iPackElem < packedToInput.size());
         packedToInput[iPackElem] = iElem;
       }
     }
   }
  }
   //get elements belonging to sampledNodes
   std::vector<int> sampleElemRanks;
   std::vector<int> sampleNodeIds(selectedNodeSet.begin(),selectedNodeSet.end()); //fill vector container with sampled nodes
   std::sort(sampleNodeIds.begin(),sampleNodeIds.end());//sort nodes in ascending order
   connections(*nodeToElem, sampleNodeIds.begin(), sampleNodeIds.end(), std::back_inserter(sampleElemRanks));//get all adjecent elements

   //fill sampled elements container
   std::vector<int> sampleElemIds;
   sampleElemIds.reserve(sampleElemRanks.size());
   std::map<int, double> weights;
   for (std::vector<int>::const_iterator it = sampleElemRanks.begin(), it_end = sampleElemRanks.end(); it != it_end; ++it) {
     const int elemRank = packedToInput[*it];
     weights.insert(std::make_pair(elemRank, 1.0));
     sampleElemIds.push_back(elemRank);
   }
 
  //construct element weight solution vector
  std::vector<double> solution(elementCount());
  for(int i = 0; i != elementCount(); i++)
    if(weights[i])
     solution[i] = weights[i];
    else 
     solution[i] = 0.0;
 
  //output Full Weights for compatibility with full mesh hyperreduction
  {
   const std::string fileName = domain->solInfo().reducedMeshFile;
   const std::ios_base::openmode mode = std::ios_base::out;
   std::ofstream weightOut(fileName.c_str(), mode);
   weightOut.precision(std::numeric_limits<double>::digits10+1);

   weightOut << "ATTRIBUTES\n";
   for (int i = 0; i != solution.size(); i++) {
    if(i == 0 ){
      if(domain->solInfo().reduceFollower) 
        weightOut << i + 1 << " 1 " << "HRC REDFOL" << " " << solution[i] << "\n";
      else
        weightOut << i + 1 << " 1 " << "HRC" << " " << solution[i] << "\n";
    } else {
      weightOut << i + 1 << " 1 " << "HRC" << " " << solution[i] << "\n";
    }
   }   
   
   weightOut << "*\n";
   weightOut << "SNSLOT\n";
   for(std::vector<std::pair<int,int> >::iterator it = compressedNodeKey.begin(); it != compressedNodeKey.end(); it++)
     weightOut << it->first + 1 << " " << it->second << std::endl; 
  }
 
  //construct and print renumbered mesh
  const FileNameInfo fileInfo;
  const MeshRenumbering meshRenumbering(sampleElemIds.begin(), sampleElemIds.end(), *elemToNode, true);
  const MeshDesc reducedMesh(domain, geoSource, meshRenumbering, weights);
  outputMeshFile(fileInfo, reducedMesh, podBasis_.vectorCount());

  // output the reduced forces
  std::ofstream meshOut(getMeshFilename(fileInfo).c_str(), std::ios_base::app);
  if(domain->solInfo().reduceFollower) meshOut << "REDFOL\n";
  meshOut << "*\nFORCES\nMODAL\n";
  meshOut.precision(std::numeric_limits<double>::digits10+1);
  for(int i=0; i<podBasis_.vectorCount(); ++i)
    meshOut << i+1 << " "  << constForceRed[i] << std::endl; 
  meshOut << "\nSNSLOT\n";
  for(std::vector<std::pair<int,int> >::iterator it = compressedNodeKey.begin(); it != compressedNodeKey.end(); it++)
   meshOut << it->first + 1 << " " << it->second << std::endl;
 
}

void
UDEIMSamplingDriver::buildForceArray(VecBasis &forceBasis, const VecBasis &displac, const VecBasis *veloc,
                                     const VecBasis *accel, std::vector<double> timeStamps_, std::vector<int> snapshotCounts_)
{//this memeber function is for converting state snapshots to force snapshots in the absence of precollected force snapshots from model I
  //most of the code is copied from assembleTrainingData in ElementSamplingDriver.C
  std::vector<double>::iterator timeStampIt = timeStamps_.begin();

  forceBasis.dimensionIs(displac.vectorCount(), displac.vectorInfo());

  FullSquareMatrix *kelArrayCopy(kelArray);

  int iSnap = 0;
  double gamma  = domain->solInfo().newmarkGamma;
  double alphaf = domain->solInfo().newmarkAlphaF;
  for(int i = 0; i < snapshotCounts_.size(); i++) {
    GenVector<double> FNLint(solVecInfo());
    GenVector<double> FLint(solVecInfo());
    for(int jSnap = 0; jSnap != snapshotCounts_[i]; ++iSnap, ++jSnap){
      filePrint(stderr,"\r %4.2f%% complete, time = %f", double(iSnap)/double(std::accumulate(snapshotCounts_.begin(),snapshotCounts_.end(),0))*100.,*timeStampIt);

      geomState->explicitUpdate(domain->getNodes(), displac[iSnap]);
      if(veloc){ geomState->setVelocity((*veloc)[iSnap]);} //just set the velocity at the nodes
      if(accel){ geomState->setAcceleration((*accel)[iSnap]);} //just set the acceleration at the nodes

      for (int iElem = 0; iElem != elementCount(); iElem++){

        //set up for internal force vector
        Vector dummy(solVecInfo()); 
        dummy = displac[iSnap];
        SingleDomainDynamic::getInternalForce( dummy, FNLint, *timeStampIt, jSnap);
        domain->getKtimesU(dummy, bcx, FLint, 1.0, kelArrayCopy);
    
        //set vector in force snapshot container
        forceBasis[iSnap] = (FNLint+FLint); 
     
        //std::cout << " NLF + LF norm = " << FNLint.norm() << " NLF norm = " << forceBasis[iSnap].norm() << " LF norm = " << FLint.norm() << std::endl;
   
        timeStampIt++;
      }
    }
  }
  filePrint(stderr,"\r %4.2f%% complete\n", 100.);
}

void UDEIMSamplingDriver::OrthoForceSnap(VecBasis &forceBasis,std::vector<double> &SVs)
{
#ifdef USE_EIGEN3
  std::cout << "... Orthogonalizting Snapshots ..." << std::endl;
  SVs.resize(forceBasis.numVectors());
  Eigen::Map<Eigen::Matrix<double,Eigen::Dynamic,1> > SingularValueMap(SVs.data(),forceBasis.numVectors());
  Eigen::Map<Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> > ForceMap(forceBasis.data(), forceBasis.size(), forceBasis.numVectors());
  Eigen::JacobiSVD<Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> > ForceSVD(ForceMap, Eigen::ComputeThinU);
  ForceMap = ForceSVD.matrixU();
  SingularValueMap = ForceSVD.singularValues();
#endif
}

int UDEIMSamplingDriver::elementCount() const {
  return domain->numElements();
}

//member functions for fascilitating the reading of snapshot data
int 
UDEIMSamplingDriver::snapSize(BasisId::Type type, std::vector<int> &snapshotCounts)
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
  return std::accumulate(snapshotCounts.begin(), snapshotCounts.end(), 0);
}

void 
UDEIMSamplingDriver::readAndProjectSnapshots(BasisId::Type type, const int vectorSize, VecBasis &podBasis,
                                              const VecNodeDof6Conversion *vecDofConversion,
                                              std::vector<int> &snapshotCounts, std::vector<double> &timeStamps, VecBasis &config)
{
  const int snapshotCount = snapSize(type, snapshotCounts);
  filePrint(stderr, " ... Reading in and Projecting %d %s Snapshots ...\n", snapshotCount, toString(type).c_str());

  config.dimensionIs(snapshotCount, vectorSize);
  timeStamps.clear();
  timeStamps.reserve(snapshotCount);

  const int skipFactor = std::max(domain->solInfo().skipPodRom, 1); // skipFactor must be >= 1
  const int skipOffSet = std::max(domain->solInfo().skipOffSet, 0); // skipOffSet must be >= 0
  const int podVectorCount = podBasis.vectorCount();
  Vector snapshot(vectorSize);
  Vector podComponents(podVectorCount);
  const FileNameInfo fileInfo;

  int offset = 0;
  for(int i = 0; i < FileNameInfo::size(type, BasisId::SNAPSHOTS); i++) {
    std::string fileName = BasisFileId(fileInfo, type, BasisId::SNAPSHOTS, i);
    filePrint(stderr, " ... Processing File: %s ...\n", fileName.c_str());
    BasisInputStream in(fileName, *vecDofConversion);

    int count = 0;
    int skipCounter = skipFactor - skipOffSet;
    while(count < snapshotCounts[i]) {
      std::pair<double, double *> data;
      data.second = snapshot.data();
      in >> data;
      assert(in);
      if(skipCounter == skipFactor) {
    //    expand(podBasis, reduce(podBasis, snapshot, podComponents), config[offset+count]);
        config[offset+count] = snapshot;
        timeStamps.push_back(data.first);
        skipCounter = 1;
        ++count;
        filePrint(stderr, "\r ... timeStamp = %7.2e, %4.2f%% complete ...", data.first, double(count)/snapshotCounts[i]*100);
      }
      else {
        ++skipCounter;
      }
    }

    filePrint(stderr,"\n");
    offset += snapshotCounts[i];
  }

  assert(timeStamps.size() == snapshotCount);
}

} /* end namespace Rom */

Rom::DriverInterface *udeimSamplingDriverNew(Domain *domain) {
  return new Rom::UDEIMSamplingDriver(domain);
}
