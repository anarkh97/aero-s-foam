#include "DistrElementSamplingDriver.h"
#include "SubElementSamplingDriver.h"
#include "RenumberingUtils.h"
#include "MeshDesc.h"
#include "DistrDomainUtils.h"

#include "DistrVecBasis.h"
#include "DistrVecBasisOps.h"

#include "DistrMasterMapping.h"
#include "DistrNodeDof6Buffer.h"
#include "DistrVecNodeDof6Conversion.h"
#include "PtrPtrIterAdapter.h"

#include "FileNameInfo.h"
#include "DistrBasisFile.h"

#include <Utils.d/DistHelper.h>

#include <algorithm>
#include <memory>
#include <cassert>
#include <iostream>

extern GeoSource *geoSource;

namespace Rom {


const DistrInfo&
DistrElementSamplingDriver::vectorSize() const {
  return decDomain->masterSolVecInfo();
}

DistrElementSamplingDriver::DistrElementSamplingDriver(Domain *domain, Communicator *comm) :
  MultiDomainDynam(domain),
  domain_(domain),
  comm_(comm)
{}

int
DistrElementSamplingDriver::snapSize(BasisId::Type type){
  std::vector<std::string> files;
  FileNameInfo fileInfo;
  if(type == BasisId::STATE) files = domain->solInfo().statePodRomFile;
  if(type == BasisId::VELOCITY) files = domain_->solInfo().velocPodRomFile;
  if(type == BasisId::ACCELERATION) files = domain_->solInfo().accelPodRomFile;
  int stateCount = 0;
  const int skipFactor = domain_->solInfo().skipPodRom;
  const int skipOffSet = domain_->solInfo().skipOffSet;
  for(int i=0; i < files.size(); i++){
    std::string fileName = BasisFileId(fileInfo,type,BasisId::SNAPSHOTS,i);
    DistrBasisInputFile in(fileName);
    stateCount += (in.stateCount() % 2) + (in.stateCount() - skipOffSet) / skipFactor;
  }
  return stateCount;
}

void
DistrElementSamplingDriver::solve() {
  //decDomain->preProcess();
  
  MultiDomainDynam::preProcess();

  FileNameInfo fileInfo;

  // Read order reduction data
  DistrVecBasis podBasis;
  DistrBasisInputFile podBasisFile(BasisFileId(fileInfo, BasisId::STATE, BasisId::POD));

  const int projectionSubspaceSize = domain_->solInfo().maxSizePodRom ?
                                     std::min(domain_->solInfo().maxSizePodRom, podBasisFile.stateCount()) :
                                     podBasisFile.stateCount();

  filePrint(stderr, " ... Projection subspace of dimension = %d ...\n", projectionSubspaceSize);
  podBasis.dimensionIs(projectionSubspaceSize, vectorSize());

  DistrVecNodeDof6Conversion converter(decDomain->getAllSubDomains(), decDomain->getAllSubDomains() + decDomain->getNumSub());

  typedef PtrPtrIterAdapter<SubDomain> SubDomIt;
  DistrMasterMapping masterMapping(SubDomIt(decDomain->getAllSubDomains()),
                                   SubDomIt(decDomain->getAllSubDomains() + decDomain->getNumSub()));
  DistrNodeDof6Buffer buffer(masterMapping.localNodeBegin(), masterMapping.localNodeEnd());
  int counter = 0;
  for (DistrVecBasis::iterator it = podBasis.begin(),
                               it_end = podBasis.end();
                               it != it_end; ++it) {
    assert(podBasisFile.validCurrentState());

    podBasisFile.currentStateBuffer(buffer);
    converter.vector(buffer, *it);

    podBasisFile.currentStateIndexInc();
    ++counter;
    filePrint(stderr,"\r %4.2f%% complete", double(counter)/double(projectionSubspaceSize)*100.);
  }
    filePrint(stderr,"\n");

  const int skipFactor = domain_->solInfo().skipPodRom;
  const int skipOffSet = domain_->solInfo().skipOffSet;
  // Read state snapshots
  DistrVecBasis snapshots;
  std::vector<double> timeStamps;
  {
    BasisId::Type type = BasisId::STATE;
    const int basisStateCount = snapSize(BasisId::STATE);
    filePrint(stderr," ... Reading in %d Displacement Snapshots ...\n",basisStateCount);

    snapshots.dimensionIs(basisStateCount, vectorSize());
    timeStamps.reserve(basisStateCount);
    int counter = 0;
    int snapshotCount = 0;
    for(int i = 0; i<domain->solInfo().statePodRomFile.size(); i++){
      std::string fileName = BasisFileId(fileInfo, BasisId::STATE, BasisId::SNAPSHOTS,i);
      std::cerr << " ...Reading in displacement snapshot: " << fileName << " ...\n";
      DistrBasisInputFile in(fileName);
      int singleBasisStateCount = (in.stateCount() % 2) + (in.stateCount() - skipOffSet) / skipFactor;
      for (DistrVecBasis::iterator it = &snapshots[snapshotCount],
          it_end = &snapshots[snapshotCount+singleBasisStateCount];
          it != it_end; ++it) {
        for(int offSet=1; offSet<=skipOffSet; ++offSet){ 
          assert(in.validCurrentState());
          in.currentStateIndexInc();}

          for(int skipCounter=1; skipCounter<=skipFactor; ++skipCounter){
            assert(in.validCurrentState());

            if(skipCounter == 1) {
              in.currentStateBuffer(buffer);
              converter.vector(buffer, *it);
              timeStamps.push_back(in.currentStateHeaderValue());
              ++counter;
              filePrint(stderr,"\rtimeStamp = %f, %4.2f%% complete",in.currentStateHeaderValue(),double(counter)/double(basisStateCount)*100.);
            }
            in.currentStateIndexInc();
          }
      }
      filePrint(stderr,"\n");
      snapshotCount += singleBasisStateCount;
    }
  }

  // Read velocity snapshots
  DistrVecBasis *velocSnapshots = 0;
  if(!domain_->solInfo().velocPodRomFile.empty()) {
    std::vector<double> timeStamps;
    velocSnapshots = new DistrVecBasis;
    const int basisStateCount = snapSize(BasisId::VELOCITY);
    filePrint(stderr," ... Reading in %d Velocity Snapshots ...\n",basisStateCount);

    velocSnapshots->dimensionIs(basisStateCount, vectorSize());
    timeStamps.reserve(basisStateCount);

    int counter = 0;
    int snapshotCount = 0;
    for(int i=0; i<domain->solInfo().velocPodRomFile.size(); i++){
      std::string fileName = BasisFileId(fileInfo, BasisId::VELOCITY, BasisId::SNAPSHOTS,i);
      std::cerr << " ...Reading in velocity snapshot: " << fileName << " ...\n";
      DistrBasisInputFile in(fileName);
      int singleBasisStateCount = (in.stateCount() % 2) + (in.stateCount() - skipOffSet) / skipFactor;
      int snapshotCount = 0;
      for (DistrVecBasis::iterator it = &((*velocSnapshots)[snapshotCount]),
          it_end = &((*velocSnapshots)[snapshotCount+singleBasisStateCount]);
          it != it_end; ++it) {
        for(int offSet=1; offSet<=skipOffSet; ++offSet){
          assert(in.validCurrentState());
          in.currentStateIndexInc();}

          for(int skipCounter=1; skipCounter<=skipFactor; ++skipCounter){
            assert(in.validCurrentState());

            if(skipCounter == 1) {
              in.currentStateBuffer(buffer);
              converter.vector(buffer, *it);
              timeStamps.push_back(in.currentStateHeaderValue());
              ++counter;
              filePrint(stderr,"\rtimeStamp = %f, %4.2f%% complete", in.currentStateHeaderValue(), double(counter)/double(basisStateCount)*100.);
            }
            in.currentStateIndexInc();
          }
      }
      filePrint(stderr,"\n");
      snapshotCount += singleBasisStateCount;
    }
  }

  // Read acceleration snapshots
  DistrVecBasis *accelSnapshots = 0;
  if(!domain_->solInfo().accelPodRomFile.empty()) {
    std::vector<double> timeStamps;
    accelSnapshots = new DistrVecBasis;
    const int basisStateCount = snapSize(BasisId::ACCELERATION);
    filePrint(stderr," ... Reading in %d Acceleration Snapshots ...\n",basisStateCount);
    accelSnapshots->dimensionIs(basisStateCount, vectorSize());
    timeStamps.reserve(basisStateCount);
  
    int skipCounter = skipFactor - skipOffSet;
    int counter = 0;
    int snapshotCount = 0;
    for(int i=0; i<domain->solInfo().accelPodRomFile.size(); i++){
      std::string fileName = BasisFileId(fileInfo, BasisId::ACCELERATION, BasisId::SNAPSHOTS,i);
      std::cerr << " ...Reading in acceleration snapshot: " << fileName << " ...\n";
      DistrBasisInputFile in(fileName);
      int singleBasisStateCount = (in.stateCount() % 2) + (in.stateCount() - skipOffSet) / skipFactor;
      for (DistrVecBasis::iterator it = &((*accelSnapshots)[snapshotCount]),
          it_end = &((*accelSnapshots)[snapshotCount+singleBasisStateCount]);
          it != it_end; ++it) {
        for(int offSet=1; offSet<=skipOffSet; ++offSet){
          assert(in.validCurrentState());
          in.currentStateIndexInc();}

          for(int skipCounter=1; skipCounter<=skipFactor; ++skipCounter){
            assert(in.validCurrentState());

            if(skipCounter == 1) {
              in.currentStateBuffer(buffer);
              converter.vector(buffer, *it);
              timeStamps.push_back(in.currentStateHeaderValue());
              ++counter;
              filePrint(stderr,"\rtimeStamp = %f, %4.2f%% complete", in.currentStateHeaderValue(),double(counter)/double(basisStateCount)*100.);
            }
            in.currentStateIndexInc();
          }
      }
      snapshotCount += singleBasisStateCount;
      filePrint(stderr,"\n");
    }
  }

  const int podVectorCount = podBasis.vectorCount();
  const int snapshotCount  = snapshots.vectorCount();

  // Temporary buffers shared by all iterations
  Vector podComponents(podVectorCount);

  //BEGIN PROJECTION-------------------------------------------------------------------------
 
  // Project snapshots on POD basis to get training configurations
  filePrint(stderr," ... Projecting displacement snapshots for training configuration ...\n");
  DistrVecBasis displac(snapshotCount, vectorSize());
  for (int iSnap = 0; iSnap != snapshotCount; ++iSnap) {
    expand(podBasis, reduce(podBasis, snapshots[iSnap], podComponents), displac[iSnap]);
    filePrint(stderr,"\r %4.2f%% complete", double(iSnap)/double(snapshotCount-1)*100.);
  }
  filePrint(stderr,"\n");

  DistrVecBasis *veloc = 0;
  if(velocSnapshots) {
    veloc = new DistrVecBasis(velocSnapshots->vectorCount(), vectorSize());

    // Project velocity snapshots on POD basis to get training configurations
    filePrint(stderr," ... Projecting velocity snapshots for training configuration ...\n");
    for (int iSnap = 0; iSnap != velocSnapshots->vectorCount(); ++iSnap) {
      expand(podBasis, reduce(podBasis, (*velocSnapshots)[iSnap], podComponents), (*veloc)[iSnap]);
      filePrint(stderr,"\r %4.2f%% complete", double(iSnap)/double(velocSnapshots->vectorCount())*100.);
    }
    filePrint(stderr,"\n");
    delete velocSnapshots;
  }

  DistrVecBasis *accel = 0;
  if(accelSnapshots) {
    accel = new DistrVecBasis(accelSnapshots->vectorCount(), vectorSize());

    // Project acceleration snapshots on POD basis to get training configurations
    filePrint(stderr," ... Projecting acceleration snapshots for training configuration ...\n");
    for (int iSnap = 0; iSnap != accelSnapshots->vectorCount(); ++iSnap) {
      expand(podBasis, reduce(podBasis, (*accelSnapshots)[iSnap], podComponents), (*accel)[iSnap]);
      filePrint(stderr,"\r %4.2f%% complete", double(iSnap)/double(accelSnapshots->vectorCount())*100.);
    }
    filePrint(stderr,"\n");
    delete accelSnapshots;
  }

  //END PROJECTION
  std::cerr << "ended projection\n";

  //Projections complete so read in normalized basis
  if(domain->solInfo().newmarkBeta == 0){
    std::string normalizedBasisFileName = BasisFileId(fileInfo,BasisId::STATE,BasisId::POD);
    normalizedBasisFileName.append(".normalized");
    DistrBasisInputFile normalizedBasisFile(normalizedBasisFileName);
    DistrVecBasis normalizedBasis;
    normalizedBasis.dimensionIs(projectionSubspaceSize,vectorSize());
    for(DistrVecBasis::iterator it = normalizedBasis.begin(), it_end = normalizedBasis.end(); it != it_end; ++it){
      assert(normalizedBasisFile.validCurrentState());
      normalizedBasisFile.currentStateBuffer(buffer);
      converter.vector(buffer, *it);
      normalizedBasisFile.currentStateIndexInc();
    }
    podBasis = normalizedBasis;
  }
  
  SubElementSamplingDriver **subDrivers = new SubElementSamplingDriver * [decDomain->getNumSub()];
  Vector *solutions = new Vector[decDomain->getNumSub()];
  int numCPUs = (structCom) ? structCom->numCPUs() : 1;
  int myID = (structCom) ? structCom->myID() : 0;
  bool verboseFlag = (myID == 0); // output to the screen only for subdomains assigned to mpi process with rank 0
#if defined(_OPENMP)
  #pragma omp parallel for schedule(static,1)
#endif
  for(int i=0; i<decDomain->getNumSub(); ++i) {
    subDrivers[i] = new SubElementSamplingDriver(decDomain->getAllSubDomains()[i]);

    VecBasis &subPodBasis = subDrivers[i]->podBasis();
    subPodBasis.dimensionIs(podVectorCount, subDrivers[i]->vectorSize());
    for(int j=0; j<podVectorCount; ++j) {
      subPodBasis[j] = StackVector(podBasis[j].subData(i), podBasis[j].subLen(i));
    }

    VecBasis &subDisplac = subDrivers[i]->displac();
    subDisplac.dimensionIs(snapshotCount, subDrivers[i]->vectorSize());
    for(int j=0; j<snapshotCount; ++j) {
      subDisplac[j] = StackVector(displac[j].subData(i), displac[j].subLen(i));
    }

    subDrivers[i]->timeStampsIs(timeStamps);
    if(veloc) {
      VecBasis *subVeloc = subDrivers[i]->veloc();
      subVeloc->dimensionIs(snapshotCount, subDrivers[i]->vectorSize());
      for(int j=0; j<snapshotCount; ++j) {
        (*subVeloc)[j] = StackVector((*veloc)[j].subData(i), (*veloc)[j].subLen(i));
      }
    }
    if(accel) {
      VecBasis *subAccel = subDrivers[i]->accel();
      subAccel->dimensionIs(snapshotCount, subDrivers[i]->vectorSize());
      for(int j=0; j<snapshotCount; ++j) {
        (*subAccel)[j] = StackVector((*accel)[j].subData(i), (*accel)[j].subLen(i));
      }
    }
    subDrivers[i]->computeSolution(solutions[i], verboseFlag);
  }
  
  std::vector<double> lweights; 
  std::vector<int> lelemIds;
  for(int i=0; i<decDomain->getNumSub(); i++) {
    subDrivers[i]->getGlobalWeights(solutions[i], lweights, lelemIds, verboseFlag);
  }
  
  std::vector<double> gweights(domain_->numElements());
  std::vector<int> gelemIds(domain_->numElements());
  
  //Gather weights and IDs from all processors
  int numLocalElems = lweights.size();
  if(structCom){
    int recvcnts[numCPUs];
    int displacements[numCPUs];
    structCom->allGather(&numLocalElems,1,&recvcnts[0],1);
    int location = 0 ;
    for(int i = 0 ; i < numCPUs ; i++){
	displacements[i]=location;
	location += recvcnts[i];
    }
    structCom->gatherv(&lweights[0], lweights.size(), &gweights[0], &recvcnts[0], &displacements[0], 0);
    structCom->gatherv(&lelemIds[0], lelemIds.size(), &gelemIds[0], &recvcnts[0], &displacements[0], 0);
  }
  else{ //no MPI
    gweights = lweights;
    gelemIds = lelemIds;
  }

/*
#ifdef USE_EIGEN3
  filePrint(stderr, " ... Compressing Basis              ...\n");
  std::vector<std::vector<int> > packedWeightedNodes(decDomain->getNumSub());
  DofSetArray **all_cdsa = new DofSetArray * [decDomain->getNumSub()];
  for(int i=0; i<decDomain->getNumSub(); ++i) {
    std::vector<int> &subWeightedNodes = packedWeightedNodes[i];

    // loop over the elements with non-zero weights and add their nodes to list
    std::vector<int> node_buffer;
    for(int j=0; j<solutions[i].size(); ++j) {
      if(solutions[i][j] > 0.0) {
        Element *ele = decDomain->getSubDomain(i)->getElementSet()[j]; 
        node_buffer.resize(ele->numNodes());
        ele->nodes(node_buffer.data());
        subWeightedNodes.insert(subWeightedNodes.end(), node_buffer.begin(), node_buffer.end());
      }
    }

    //sort nodes in ascending order and erase redundant nodes
    std::sort(subWeightedNodes.begin(), subWeightedNodes.end());
    std::vector<int>::iterator packedNodeIt = std::unique(subWeightedNodes.begin(),subWeightedNodes.end());
    subWeightedNodes.resize(packedNodeIt-subWeightedNodes.begin());

    all_cdsa[i] = decDomain->getSubDomain(i)->getCDSA();
  }
  podBasis.makeSparseBasis(packedWeightedNodes, all_cdsa);
  delete [] all_cdsa;
  // TODO print compressed basis to file
#endif
*/
  if(myID==0){
     //Weights output file generation
     const std::string fileName = domain_->solInfo().reducedMeshFile;
     std::ofstream weightOut(fileName.c_str(),std::ios_base::out);
     weightOut << "ATTRIBUTES\n";
     for(int i = 0 ; i < gweights.size(); i++){
	weightOut<< gelemIds[i]+1 << " 1 " << "HRC" << " " << gweights[i] << "\n";
     }  

    //Mesh output file generation
    std::map<int,double> weightsMap;
    std::vector<int> reducedelemIds;
    for(int i =0 ; i< gweights.size(); i++){
      if(gweights[i]>0){
	reducedelemIds.push_back(gelemIds[i]);
        weightsMap.insert(std::pair<int,double>(gelemIds[i],gweights[i]));
      }
    }

    const FileNameInfo fileInfo;
    for(int i=0; i<decDomain->getNumSub(); i++) {
      decDomain->getSubDomain(i)->renumberElementsGlobal();
    }

    Elemset &inputElemSet = *(geoSource->getElemSet());
    std::auto_ptr<Connectivity> elemToNode(new Connectivity(&inputElemSet));

    const MeshRenumbering meshRenumbering(reducedelemIds.begin(), reducedelemIds.end(), *elemToNode, verboseFlag);
    const MeshDesc reducedMesh(domain_, geoSource, meshRenumbering, weightsMap); 
    outputMeshFile(fileInfo, reducedMesh, podBasis.vectorCount());
  }

  if(structCom) structCom->sync();

  delete [] subDrivers;
  delete [] solutions;
  if(veloc) delete veloc;
  if(accel) delete accel;
}

} /* end namespace Rom */

extern Communicator *structCom;

Rom::DriverInterface *distrElementSamplingDriverNew(Domain *domain) {
  return new Rom::DistrElementSamplingDriver(domain, structCom);
}
