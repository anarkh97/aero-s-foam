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

#include "BasisFileStream.h"
#include "VecBasisFile.h"

#include <Utils.d/DistHelper.h>

#include <algorithm>
#include <memory>
#include <cassert>
#include <iostream>

extern GeoSource *geoSource;

namespace Rom {

std::string getMeshFilename(const FileNameInfo &fileInfo);

const DistrInfo&
DistrElementSamplingDriver::vectorSize() const {
  return decDomain->masterSolVecInfo();
}

DistrElementSamplingDriver::DistrElementSamplingDriver(Domain *domain, Communicator *comm) :
  MultiDomainDynam(domain),
  comm_(comm)
{}

int
DistrElementSamplingDriver::snapSize(BasisId::Type type) {
  std::vector<std::string> files;
  FileNameInfo fileInfo;
  if(type == BasisId::STATE) files = domain->solInfo().statePodRomFile;
  if(type == BasisId::VELOCITY) files = domain->solInfo().velocPodRomFile;
  if(type == BasisId::ACCELERATION) files = domain->solInfo().accelPodRomFile;
  int stateCount = 0;
  const int skipFactor = domain->solInfo().skipPodRom;
  const int skipOffSet = domain->solInfo().skipOffSet;
  for(int i=0; i < files.size(); i++){
    std::string fileName = BasisFileId(fileInfo,type,BasisId::SNAPSHOTS,i);
    DistrBasisInputFile in(fileName);
    stateCount += (in.stateCount() % 2) + (in.stateCount() - skipOffSet) / skipFactor;
  }
  return stateCount;
}

void
DistrElementSamplingDriver::solve() {
  
  MultiDomainDynam::preProcess();

  FileNameInfo fileInfo;

  // Read order reduction data
  DistrVecBasis podBasis;
  DistrBasisInputFile podBasisFile(BasisFileId(fileInfo, BasisId::STATE, BasisId::POD));

  const int projectionSubspaceSize = domain->solInfo().maxSizePodRom ?
                                     std::min(domain->solInfo().maxSizePodRom, podBasisFile.stateCount()) :
                                     podBasisFile.stateCount();

  //filePrint(stderr, " ... Projection subspace of dimension = %d ...\n", projectionSubspaceSize);
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

  const int skipFactor = domain->solInfo().skipPodRom;
  const int skipOffSet = domain->solInfo().skipOffSet;
  // Read state snapshots
  DistrVecBasis snapshots;
  std::vector<double> timeStamps;
  {
    BasisId::Type type = BasisId::STATE;
    const int basisStateCount = snapSize(BasisId::STATE);
    filePrint(stderr, " ... Reading in %d Displacement Snapshots ...\n", basisStateCount);

    snapshots.dimensionIs(basisStateCount, vectorSize());
    timeStamps.reserve(basisStateCount);

    int counter = 0;
    int snapshotCount = 0;
    for(int i = 0; i < domain->solInfo().statePodRomFile.size(); i++) {
      std::string fileName = BasisFileId(fileInfo, BasisId::STATE, BasisId::SNAPSHOTS,i);
      filePrint(stderr, " ... Processing File: %s ...\n", fileName.c_str());
      DistrBasisInputFile in(fileName);
      int singleBasisStateCount = (in.stateCount() % 2) + (in.stateCount() - skipOffSet) / skipFactor;
      for(DistrVecBasis::iterator it = &snapshots[snapshotCount],
          it_end = &snapshots[snapshotCount+singleBasisStateCount];
          it != it_end; ++it) {
        for(int offSet = 1; offSet <= skipOffSet; ++offSet) {
          assert(in.validCurrentState());
          in.currentStateIndexInc();
        }

        for(int skipCounter = 1; skipCounter <= skipFactor; ++skipCounter) {
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
  if(!domain->solInfo().velocPodRomFile.empty()) {
    std::vector<double> timeStamps;
    velocSnapshots = new DistrVecBasis;
    const int basisStateCount = snapSize(BasisId::VELOCITY);
    filePrint(stderr, " ... Reading in %d Velocity Snapshots ...\n", basisStateCount);

    velocSnapshots->dimensionIs(basisStateCount, vectorSize());
    timeStamps.reserve(basisStateCount);

    int counter = 0;
    int snapshotCount = 0;
    for(int i = 0; i < domain->solInfo().velocPodRomFile.size(); i++) {
      std::string fileName = BasisFileId(fileInfo, BasisId::VELOCITY, BasisId::SNAPSHOTS,i);
      filePrint(stderr, " ... Processing File: %s ...\n", fileName.c_str());
      DistrBasisInputFile in(fileName);
      int singleBasisStateCount = (in.stateCount() % 2) + (in.stateCount() - skipOffSet) / skipFactor;
      for(DistrVecBasis::iterator it = &((*velocSnapshots)[snapshotCount]),
          it_end = &((*velocSnapshots)[snapshotCount+singleBasisStateCount]);
          it != it_end; ++it) {
        for(int offSet = 1; offSet <= skipOffSet; ++offSet) {
          assert(in.validCurrentState());
          in.currentStateIndexInc();
        }

        for(int skipCounter = 1; skipCounter <= skipFactor; ++skipCounter) {
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

  // Read acceleration snapshots
  DistrVecBasis *accelSnapshots = 0;
  if(!domain->solInfo().accelPodRomFile.empty()) {
    std::vector<double> timeStamps;
    accelSnapshots = new DistrVecBasis;
    const int basisStateCount = snapSize(BasisId::ACCELERATION);
    filePrint(stderr, " ... Reading in %d Acceleration Snapshots ...\n", basisStateCount);

    accelSnapshots->dimensionIs(basisStateCount, vectorSize());
    timeStamps.reserve(basisStateCount);

    int counter = 0;
    int snapshotCount = 0;
    for(int i = 0; i < domain->solInfo().accelPodRomFile.size(); i++){
      std::string fileName = BasisFileId(fileInfo, BasisId::ACCELERATION, BasisId::SNAPSHOTS, i);
      filePrint(stderr, " ... Processing File: %s ...\n", fileName.c_str());
      DistrBasisInputFile in(fileName);
      int singleBasisStateCount = (in.stateCount() % 2) + (in.stateCount() - skipOffSet) / skipFactor;
      for(DistrVecBasis::iterator it = &((*accelSnapshots)[snapshotCount]),
          it_end = &((*accelSnapshots)[snapshotCount+singleBasisStateCount]);
          it != it_end; ++it) {
        for(int offSet = 1; offSet <= skipOffSet; ++offSet) {
          assert(in.validCurrentState());
          in.currentStateIndexInc();
        }

        for(int skipCounter = 1; skipCounter <= skipFactor; ++skipCounter) {
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
      filePrint(stderr,"\r %4.2f%% complete", double(iSnap)/double(velocSnapshots->vectorCount()-1)*100.);
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
      filePrint(stderr,"\r %4.2f%% complete", double(iSnap)/double(accelSnapshots->vectorCount()-1)*100.);
    }
    filePrint(stderr,"\n");
    delete accelSnapshots;
  }

  //END PROJECTION

  //Projections complete so read in normalized basis
  if(domain->solInfo().newmarkBeta == 0) {
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
 
  int glNumSubs = decDomain->getNumSub();
  structCom->globalSum(1,&glNumSubs);
 
  SubElementSamplingDriver **subDrivers = new SubElementSamplingDriver * [decDomain->getNumSub()];
  Vector *solutions = new Vector[decDomain->getNumSub()];
  Vector *trainingTargets = new Vector[decDomain->getNumSub()];
  double *targetMagnitudes = new double[decDomain->getNumSub()];
  int numCPUs = (structCom) ? structCom->numCPUs() : 1;
  int myID = (structCom) ? structCom->myID() : 0;
  bool verboseFlag = (myID == 0); // output to the screen only for subdomains assigned to mpi process with rank 0
  double glTargMagnitude = 0;
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

    subDrivers[i]->preProcess();

    trainingTargets[i].reset(subPodBasis.vectorCount()*subDisplac.vectorCount(), 0.0);
    subDrivers[i]->assembleTrainingData(trainingTargets[i]);
    targetMagnitudes[i] = norm(trainingTargets[i]);
    glTargMagnitude += targetMagnitudes[i]*targetMagnitudes[i];
  }

  if(structCom)
    structCom->globalSum(1,&glTargMagnitude);
  glTargMagnitude = sqrt(glTargMagnitude);

#if defined(_OPENMP)
  #pragma omp parallel for schedule(static,1)
#endif
  for(int i=0; i<decDomain->getNumSub(); ++i) {
    double relativeTolerance = (domain->solInfo().localTol) ? domain->solInfo().tolPodRom*glTargMagnitude/(glNumSubs*targetMagnitudes[i])
                                                            : domain->solInfo().tolPodRom;
    if(verboseFlag) filePrint(stderr, " ... Training Tolerance for SubDomain %d is %f ...\n", decDomain->getSubDomain(i)->subNum()+1, relativeTolerance);
    subDrivers[i]->computeSolution(trainingTargets[i], solutions[i], relativeTolerance, verboseFlag);
  }
  delete [] trainingTargets;
  delete [] targetMagnitudes;
  
  std::vector<double> lweights; 
  std::vector<int> lelemIds;
  for(int i=0; i<decDomain->getNumSub(); i++) {
    subDrivers[i]->getGlobalWeights(solutions[i], lweights, lelemIds, verboseFlag);
  }
  
  std::vector<double> gweights(domain->numElements());
  std::vector<int> gelemIds(domain->numElements());
  
  //Gather weights and IDs from all processors
  int numLocalElems = lweights.size();
  if(structCom) {
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

  // compute the reduced forces (constant only)
  DistrVector constForceFull(MultiDomainDynam::solVecInfo());
  MultiDomainDynam::getConstForce(constForceFull);
  Vector constForceRed(podBasis.vectorCount());
  reduce(podBasis, constForceFull, constForceRed);

  if(myID == 0) {
     //Weights output file generation
     const std::string fileName = domain->solInfo().reducedMeshFile;
     std::ofstream weightOut(fileName.c_str(), std::ios_base::out);
     weightOut << "ATTRIBUTES\n";
     bool firstTime = true;
     for(int i = 0 ; i < gweights.size(); i++) {
       if(domain->solInfo().reduceFollower && firstTime) {
         weightOut<< gelemIds[i]+1 << " 1 " << "HRC REDFOL" << " " << gweights[i] << "\n";
         firstTime = false;
       }
       else {
	 weightOut<< gelemIds[i]+1 << " 1 " << "HRC" << " " << gweights[i] << "\n";
       }
     }  

    //Mesh output file generation
    std::map<int,double> weightsMap;
    std::vector<int> reducedelemIds;
    for(int i = 0; i < gweights.size(); i++) {
      if(gweights[i] > 0) {
	reducedelemIds.push_back(gelemIds[i]);
        weightsMap.insert(std::pair<int,double>(gelemIds[i], gweights[i]));
      }
    }

    const FileNameInfo fileInfo;
    for(int i = 0; i < decDomain->getNumSub(); i++) {
      decDomain->getSubDomain(i)->renumberElementsGlobal();
    }

    // read in the truncated basis into a (non-distributed) VecBasis
    domain->preProcessing();
    buildDomainCdsa();
    std::string fileName2 = BasisFileId(fileInfo, BasisId::STATE, BasisId::POD);
    if(domain->solInfo().newmarkBeta == 0) fileName2.append(".normalized");
    const VecNodeDof6Conversion vecDofConversion(*domain->getCDSA());
    BasisInputStream in(fileName2, vecDofConversion) ;
    VecBasis podBasis;
    const int podSizeMax = domain->solInfo().maxSizePodRom;
    if(podSizeMax != 0) {
      readVectors(in, podBasis, podSizeMax);
    } else {
      readVectors(in, podBasis);
    }

    // output the reduced mesh
    Elemset &inputElemSet = *(geoSource->getElemSet());
    std::auto_ptr<Connectivity> elemToNode(new Connectivity(&inputElemSet));
    const MeshRenumbering meshRenumbering(reducedelemIds.begin(), reducedelemIds.end(), *elemToNode, verboseFlag);
    const MeshDesc reducedMesh(domain, geoSource, meshRenumbering, weightsMap); 
    outputMeshFile(fileInfo, reducedMesh, podBasis.vectorCount());

    // output the reduced forces
    std::ofstream meshOut(getMeshFilename(fileInfo).c_str(), std::ios_base::app);
    if(domain->solInfo().reduceFollower) meshOut << "REDFOL\n";
    meshOut << "*\nFORCES\nMODAL\n";
    meshOut.precision(std::numeric_limits<double>::digits10+1);
    for(int i=0; i<podBasis.vectorCount(); ++i)
      meshOut << i+1 << " "  << constForceRed[i] << std::endl;

#ifdef USE_EIGEN3
    // build and output compressed basis
    podBasis.makeSparseBasis(meshRenumbering.reducedNodeIds(), domain->getCDSA());
    {
      std::string filename = BasisFileId(fileInfo, BasisId::STATE, BasisId::POD);
      filename.append(".reduced");
      if(domain->solInfo().newmarkBeta == 0) filename.append(".normalized");
      filePrint(stderr," ... Writing compressed basis to file %s ...\n", filename.c_str());
      DofSetArray reduced_dsa(reducedMesh.nodes().size(), const_cast<Elemset&>(reducedMesh.elements()));
      int num_bc = reducedMesh.dirichletBConds().size();
      BCond *bc = (num_bc > 0) ? const_cast<BCond*>(&reducedMesh.dirichletBConds()[0]) : NULL;
      ConstrainedDSA reduced_cdsa(reduced_dsa, num_bc, bc);
      VecNodeDof6Conversion converter(reduced_cdsa);
      BasisOutputStream output(filename, converter, false);

      for (int iVec = 0; iVec < podBasis.vectorCount(); ++iVec) {
        output << podBasis.getCompressedBasis().col(iVec);
      }
    }
#endif
  }

  if(structCom) structCom->sync();

  delete [] subDrivers;
  delete [] solutions;
  if(veloc) delete veloc;
  if(accel) delete accel;
}

void
DistrElementSamplingDriver::buildDomainCdsa() {
  const int numdof = domain->numdof();
  SimpleBuffer<int> bc(numdof);
  SimpleBuffer<double> bcx(numdof);

  domain->make_bc(bc.array(), bcx.array());
  domain->make_constrainedDSA(bc.array());
}

} /* end namespace Rom */

extern Communicator *structCom;

Rom::DriverInterface *distrElementSamplingDriverNew(Domain *domain) {
  return new Rom::DistrElementSamplingDriver(domain, structCom);
}
