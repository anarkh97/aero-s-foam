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

  // Read state snapshots
  DistrVecBasis snapshots;
  std::vector<double> timeStamps;
  {
    DistrBasisInputFile in(BasisFileId(fileInfo, BasisId::STATE, BasisId::SNAPSHOTS));
    const int skipFactor = domain_->solInfo().skipPodRom;
    const int skipOffSet = domain_->solInfo().skipOffSet;
    const int basisStateCount = (in.stateCount() % 2) + (in.stateCount() - skipOffSet) / skipFactor;
    filePrint(stderr," ... Reading in %d Displacement Snapshots ...\n",basisStateCount);

    snapshots.dimensionIs(basisStateCount, vectorSize());
    timeStamps.reserve(basisStateCount);

    int counter = 0;
    for (DistrVecBasis::iterator it = snapshots.begin(),
                                 it_end = snapshots.end();
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
  }

  // Read velocity snapshots
  DistrVecBasis *velocSnapshots = 0;
  if(domain_->solInfo().velocPodRomFile != "") {
    std::vector<double> timeStamps;
    velocSnapshots = new DistrVecBasis;
    DistrBasisInputFile in(BasisFileId(fileInfo, BasisId::VELOCITY, BasisId::SNAPSHOTS));
    const int skipFactor = domain_->solInfo().skipPodRom;
    const int skipOffSet = domain_->solInfo().skipOffSet;
    const int basisStateCount = (in.stateCount() % 2) + (in.stateCount() - skipOffSet) / skipFactor;
    filePrint(stderr," ... Reading in %d Velocity Snapshots ...\n",basisStateCount);

    velocSnapshots->dimensionIs(basisStateCount, vectorSize());
    timeStamps.reserve(basisStateCount);

    int counter = 0;
    for (DistrVecBasis::iterator it = velocSnapshots->begin(),
                                 it_end = velocSnapshots->end();
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
          filePrint(stderr,"\r %4.2f%% complete", double(counter)/double(basisStateCount)*100.);
          }
        in.currentStateIndexInc();
      }
    }
    filePrint(stderr,"\n");
  }

  // Read acceleration snapshots
  DistrVecBasis *accelSnapshots = 0;
  if(domain_->solInfo().accelPodRomFile != "") {
    std::vector<double> timeStamps;
    accelSnapshots = new DistrVecBasis;
    DistrBasisInputFile in(BasisFileId(fileInfo, BasisId::ACCELERATION, BasisId::SNAPSHOTS));
    const int skipFactor = domain_->solInfo().skipPodRom;
    const int skipOffSet = domain_->solInfo().skipOffSet;
    const int basisStateCount = (in.stateCount() % 2) + (in.stateCount() - skipOffSet) / skipFactor;
    filePrint(stderr," ... Reading in %d Acceleration Snapshots ...\n",basisStateCount);
 
    accelSnapshots->dimensionIs(basisStateCount, vectorSize());
    timeStamps.reserve(basisStateCount);
  
    int skipCounter = skipFactor - skipOffSet;
    int counter = 0;
    for (DistrVecBasis::iterator it = accelSnapshots->begin(),
                                 it_end = accelSnapshots->end();
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
          filePrint(stderr,"\r %4.2f%% complete", double(counter)/double(basisStateCount)*100.);
          }
        in.currentStateIndexInc();
      }
    }
    filePrint(stderr,"\n");
    filePrint(stderr,"counter = %d\n",counter);
  }

  const int podVectorCount = podBasis.vectorCount();
  const int snapshotCount  = snapshots.vectorCount();

  // Temporary buffers shared by all iterations
  Vector podComponents(podVectorCount);

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
    subDrivers[i]->getGlobalWeights(solutions[i], lweights, lelemIds, true, verboseFlag);
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
    outputMeshFile(fileInfo, reducedMesh);
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
