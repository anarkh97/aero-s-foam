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

// Forward declarations
// ====================
int snapSize(BasisId::Type type, std::vector<int> &snapshotCounts);

// Non-member functions
// ====================
std::string getMeshFilename(const FileNameInfo &fileInfo);

void readAndProjectSnapshots(BasisId::Type type, const DistrInfo &vectorSize, DistrNodeDof6Buffer &buffer,
                             DistrVecNodeDof6Conversion &converter, DistrVecBasis &podBasis,
                             std::vector<int> &snapshotCounts, std::vector<double> &timeStamps, DistrVecBasis *&config)
{
  const int snapshotCount = snapSize(type, snapshotCounts);
  filePrint(stderr, " ... Reading in and Projecting %d %s Snapshots ...\n", snapshotCount, toString(type).c_str());

  config = new DistrVecBasis(snapshotCount, vectorSize);
  timeStamps.clear();
  timeStamps.reserve(snapshotCount);

  const int skipFactor = std::max(domain->solInfo().skipPodRom, 1); // skipFactor must be >= 1
  const int skipOffSet = std::max(domain->solInfo().skipOffSet, 0); // skipOffSet must be >= 0
  const int podVectorCount = podBasis.vectorCount();
  DistrVector snapshot(vectorSize);
  Vector podComponents(podVectorCount);
  FileNameInfo fileInfo;

  int offset = 0;
  for(int i = 0; i < FileNameInfo::size(type, BasisId::SNAPSHOTS); i++) {
    std::string fileName = BasisFileId(fileInfo, type, BasisId::SNAPSHOTS, i);
    filePrint(stderr, " ... Processing File: %s ...\n", fileName.c_str());
    DistrBasisInputFile in(fileName);

    int count = 0;
    int skipCounter = skipFactor - skipOffSet;
    while(count < snapshotCounts[i]) {
      if(skipCounter == skipFactor) {
        assert(in.validCurrentState());
        in.currentStateBuffer(buffer);
        converter.vector(buffer, snapshot);
        expand(podBasis, reduce(podBasis, snapshot, podComponents), (*config)[offset+count]); // do projection
        timeStamps.push_back(in.currentStateHeaderValue());
        skipCounter = 1;
        ++count;
        filePrint(stderr, "\r ... timeStamp = %7.2e, %4.2f%% complete ...", in.currentStateHeaderValue(), double(count)/snapshotCounts[i]*100);
      }
      else {
        ++skipCounter;
      }
      in.currentStateIndexInc();
    }

    filePrint(stderr,"\n");
    offset += snapshotCounts[i];
  }

  assert(timeStamps.size() == snapshotCount);
}

// Member functions
// ================
const DistrInfo&
DistrElementSamplingDriver::vectorSize() const
{
  return decDomain->masterSolVecInfo();
}

DistrElementSamplingDriver::DistrElementSamplingDriver(Domain *domain, Communicator *comm) :
  MultiDomainDynam(domain),
  comm_(comm)
{}

void
DistrElementSamplingDriver::solve()
{
  MultiDomainDynam::preProcess();
  if(domain->solInfo().newmarkBeta == 0) {
    decDomain->assembleNodalInertiaTensors(melArray);
  }

  FileNameInfo fileInfo;

  DistrVecBasis podBasis;
  DistrBasisInputFile podBasisFile(BasisFileId(fileInfo, BasisId::STATE, BasisId::POD));

  const int projectionSubspaceSize = domain->solInfo().maxSizePodRom ?
                                     std::min(domain->solInfo().maxSizePodRom, podBasisFile.stateCount()) :
                                     podBasisFile.stateCount();

  filePrint(stderr, " ... Proj. Subspace Dimension = %-3d ...\n", projectionSubspaceSize);
  podBasis.dimensionIs(projectionSubspaceSize, vectorSize());

  DistrVecNodeDof6Conversion converter(decDomain->getAllSubDomains(), decDomain->getAllSubDomains() + decDomain->getNumSub());

  typedef PtrPtrIterAdapter<SubDomain> SubDomIt;
  DistrMasterMapping masterMapping(SubDomIt(decDomain->getAllSubDomains()),
                                   SubDomIt(decDomain->getAllSubDomains() + decDomain->getNumSub()));
  DistrNodeDof6Buffer buffer(masterMapping.localNodeBegin(), masterMapping.localNodeEnd());

  for (DistrVecBasis::iterator it = podBasis.begin(),
                               it_end = podBasis.end();
                               it != it_end; ++it) {
    assert(podBasisFile.validCurrentState());

    podBasisFile.currentStateBuffer(buffer);
    converter.vector(buffer, *it);

    podBasisFile.currentStateIndexInc();
  }

  const int podVectorCount = podBasis.vectorCount();

  // Read some displacement snapshots from one or more files and project them on to the basis
  DistrVecBasis *displac; 
  std::vector<double> timeStamps;
  std::vector<int> snapshotCounts;
  readAndProjectSnapshots(BasisId::STATE, vectorSize(), buffer, converter, podBasis,
                          snapshotCounts, timeStamps, displac);

  // Optionally, read some velocity snapshots and project them on to the reduced order basis
  DistrVecBasis *veloc = 0;
  if(!domain->solInfo().velocPodRomFile.empty()) {
    std::vector<double> velTimeStamps;
    std::vector<int> velSnapshotCounts;
    readAndProjectSnapshots(BasisId::VELOCITY, vectorSize(), buffer, converter, podBasis,
                            velSnapshotCounts, velTimeStamps, veloc);
    if(velSnapshotCounts != snapshotCounts) std::cerr << " *** WARNING: inconsistent velocity snapshots\n";
  }

  // Optionally, read some acceleration snapshots and project them on to the reduced order basis
  DistrVecBasis *accel = 0;
  if(!domain->solInfo().accelPodRomFile.empty()) {
    std::vector<double> accTimeStamps;
    std::vector<int> accSnapshotCounts;
    readAndProjectSnapshots(BasisId::ACCELERATION, vectorSize(), buffer, converter, podBasis,
                            accSnapshotCounts, accTimeStamps, accel);
    if(accSnapshotCounts != snapshotCounts) std::cerr << " *** WARNING: inconsistent acceleration snapshots\n";
  }

  const int snapshotCount = std::accumulate(snapshotCounts.begin(), snapshotCounts.end(), 0);

  // read in mass-normalized basis
  if(domain->solInfo().newmarkBeta == 0 || domain->solInfo().useMassNormalizedBasis) {
    std::string normalizedBasisFileName = BasisFileId(fileInfo, BasisId::STATE, BasisId::POD);
    normalizedBasisFileName.append(".normalized");
    DistrBasisInputFile normalizedBasisFile(normalizedBasisFileName);
    for(DistrVecBasis::iterator it = podBasis.begin(), it_end = podBasis.end(); it != it_end; ++it) {
      assert(normalizedBasisFile.validCurrentState());
      normalizedBasisFile.currentStateBuffer(buffer);
      converter.vector(buffer, *it);
      normalizedBasisFile.currentStateIndexInc();
    }
  }
 
  int glNumSubs = decDomain->getNumSub();
  structCom->globalSum(1,&glNumSubs);
 
  SubElementSamplingDriver **subDrivers = new SubElementSamplingDriver * [decDomain->getNumSub()];
  Vector *solutions = new Vector[decDomain->getNumSub()];
  double *targetMagnitudes = new double[decDomain->getNumSub()];
  int numCPUs = (structCom) ? structCom->numCPUs() : 1;
  int myID = (structCom) ? structCom->myID() : 0;
  bool verboseFlag = (myID == 0); // output to the screen only for subdomains assigned to mpi process with rank 0
  double glTargMagnitude = 0;
  Vector glTrainingTarget(podVectorCount*snapshotCount, 0.0);
#if defined(_OPENMP)
  #pragma omp parallel for schedule(static,1)
#endif
  for(int i = 0; i < decDomain->getNumSub(); ++i) {
    subDrivers[i] = new SubElementSamplingDriver(decDomain->getAllSubDomains()[i]);

    std::vector<StackVector> subPodBasis(podVectorCount);
    for(int j = 0; j < podVectorCount; ++j) {
      subPodBasis[j].setData(podBasis[j].subData(i), podBasis[j].subLen(i));
    }

    std::vector<StackVector> subDisplac(snapshotCount);
    for(int j = 0; j < snapshotCount; ++j) {
      subDisplac[j].setData((*displac)[j].subData(i), (*displac)[j].subLen(i));
    }

    subDrivers[i]->timeStampsIs(timeStamps);
    subDrivers[i]->snapshotCountsIs(snapshotCounts);

    std::vector<StackVector> *subVeloc = 0;
    if(veloc) {
      subVeloc = new std::vector<StackVector>(snapshotCount);
      for(int j = 0; j < snapshotCount; ++j) {
        (*subVeloc)[j].setData((*veloc)[j].subData(i), (*veloc)[j].subLen(i));
      }
    }

    std::vector<StackVector> *subAccel = 0;
    if(accel) {
      subAccel = new std::vector<StackVector>(snapshotCount); 
      for(int j = 0; j < snapshotCount; ++j) {
        (*subAccel)[j].setData((*accel)[j].subData(i), (*accel)[j].subLen(i));
      }
    }

    subDrivers[i]->preProcess();
    subDrivers[i]->solver().problemSizeIs(podVectorCount*snapshotCount, subDrivers[i]->elementCount());

    for(int j=0; j<subDrivers[i]->solver().equationCount(); ++j) subDrivers[i]->solver().rhsBuffer()[j] = 0.0;
    subDrivers[i]->assembleTrainingData(subPodBasis, podVectorCount, subDisplac, subVeloc, subAccel);
    if(!domain->solInfo().globalErr) subDrivers[i]->clean();
    if(subVeloc) delete subVeloc;
    if(subAccel) delete subAccel;

    StackVector trainingTarget(subDrivers[i]->solver().rhsBuffer(), podVectorCount*snapshotCount);
    targetMagnitudes[i] = trainingTarget.norm();
#if defined(_OPENMP)
    #pragma omp critical
#endif
    glTrainingTarget += trainingTarget;
  }
  if(!domain->solInfo().globalErr) {
    delete displac;
    if(veloc) delete veloc;
    if(accel) delete accel;
  }

  if(structCom) 
    structCom->globalSum(glTrainingTarget.size(), glTrainingTarget.data()); 
  glTargMagnitude = glTrainingTarget.norm();
  //if(structCom->myID() == 0) std::cerr << "glTargMagnitude = " << glTargMagnitude << std::endl;

  Vector glAx(podVectorCount*snapshotCount, 0.0);

#if defined(_OPENMP)
  #pragma omp parallel for schedule(static,1)
#endif
  for(int i = 0; i < decDomain->getNumSub(); ++i) {
    double relativeTolerance = (domain->solInfo().localTol) 
                               ? domain->solInfo().tolPodRom*glTargMagnitude/(glNumSubs*targetMagnitudes[i])
                               : domain->solInfo().tolPodRom;
    if(verboseFlag) filePrint(stderr, " ... Training Tolerance for SubDomain %d is %f ...\n",
                              decDomain->getSubDomain(i)->subNum()+1, relativeTolerance);
    subDrivers[i]->computeSolution(solutions[i], relativeTolerance, verboseFlag);

    if(domain->solInfo().globalErr) {
      std::vector<StackVector> subPodBasis(podVectorCount);
      for(int j = 0; j < podVectorCount; ++j) {
        subPodBasis[j].setData(podBasis[j].subData(i), podBasis[j].subLen(i));
      }

      std::vector<StackVector> subDisplac(snapshotCount);
      for(int j = 0; j < snapshotCount; ++j) {
        subDisplac[j].setData((*displac)[j].subData(i), (*displac)[j].subLen(i));
      }

      subDrivers[i]->timeStampsIs(timeStamps);
      subDrivers[i]->snapshotCountsIs(snapshotCounts);

      std::vector<StackVector> *subVeloc = 0;
      if(veloc) {
        subVeloc = new std::vector<StackVector>(snapshotCount);
        for(int j = 0; j < snapshotCount; ++j) {
          (*subVeloc)[j].setData((*veloc)[j].subData(i), (*veloc)[j].subLen(i));
        }
      }

      std::vector<StackVector> *subAccel = 0;
      if(accel) {
        subAccel = new std::vector<StackVector>(snapshotCount);
        for(int j = 0; j < snapshotCount; ++j) {
          (*subAccel)[j].setData((*accel)[j].subData(i), (*accel)[j].subLen(i));
        }
      }

      StackVector Ax(subDrivers[i]->solver().rhsBuffer(), podVectorCount*snapshotCount);
      Ax.zero();
      subDrivers[i]->assembleWeightedTrainingData(subPodBasis, podVectorCount, subDisplac, subVeloc, subAccel, solutions[i]);
      if(verboseFlag) std::cerr << " ... Norm of training solution (A*x) for subdomain " << decDomain->getSubDomain(i)->subNum()+1 
                                << " is " << Ax.norm() << " ...\n";
      subDrivers[i]->clean();
      if(subVeloc) delete subVeloc;
      if(subAccel) delete subAccel;

#if defined(_OPENMP)
      #pragma omp critical
#endif
      glAx += Ax;
    }
  }
  if(domain->solInfo().globalErr) {
    if(structCom)
      structCom->globalSum(glAx.size(), glAx.data());
    glAx -= glTrainingTarget;
    if(structCom->myID() == 0) std::cerr << " ... Global relative residual = " << glAx.norm()/glTargMagnitude << " ...\n";

    delete displac;
    if(veloc) delete veloc;
    if(accel) delete accel;
  }
  delete [] targetMagnitudes;
  
  std::vector<double> lweights; 
  std::vector<int> lelemIds;
  for(int i = 0; i < decDomain->getNumSub(); i++) {
    subDrivers[i]->getGlobalWeights(solutions[i], lweights, lelemIds, verboseFlag);
    delete subDrivers[i];
  }
  delete [] solutions;
  delete [] subDrivers;
  
  std::vector<double> gweights(domain->numElements());
  std::vector<int> gelemIds(domain->numElements());
  
  // Gather weights and IDs from all processors
  int numLocalElems = lweights.size();
  if(structCom) {
    int recvcnts[numCPUs];
    int displacements[numCPUs];
    structCom->allGather(&numLocalElems, 1, &recvcnts[0], 1);
    int location = 0;
    for(int i = 0; i < numCPUs; i++) {
      displacements[i] = location;
      location += recvcnts[i];
    }
    structCom->gatherv(&lweights[0], lweights.size(), &gweights[0], &recvcnts[0], &displacements[0], 0);
    structCom->gatherv(&lelemIds[0], lelemIds.size(), &gelemIds[0], &recvcnts[0], &displacements[0], 0);
  }
  else {
    gweights = lweights;
    gelemIds = lelemIds;
  }

  // Compute the reduced forces
  DistrVector forceFull(decDomain->masterSolVecInfo());
  GenAssembler<double> * assembler = decDomain->getSolVecAssembler();
  // 1) gravity
  Vector gravForceRed(podBasis.vectorCount());
  if(domain->gravityFlag()) {
    MultiDomainDynam::getGravityForce(forceFull);
    assembler->assemble(forceFull);
    reduce(podBasis, forceFull, gravForceRed);
  }
  // 2) constant force or constant part of time-dependent forces (default loadset only) TODO add support for multiple loadsets
  MultiDomainDynam::getUnamplifiedExtForce(forceFull, 0);
  assembler->assemble(forceFull);
  Vector constForceRed(podBasis.vectorCount());
  bool reduce_f = (forceFull.norm() != 0);
  if(reduce_f) reduce(podBasis, forceFull, constForceRed);

  // compute the reduced initial conditions
  DistrVector d0Full(MultiDomainDynam::solVecInfo()),
              v0Full(MultiDomainDynam::solVecInfo());
  DistrVector tmp(decDomain->masterSolVecInfo());
  SysState<DistrVector> inState(d0Full, v0Full, tmp, tmp);
  MultiDomainDynam::getInitState(inState);
  Vector d0Red(podBasis.vectorCount()),
         v0Red(podBasis.vectorCount());
  bool reduce_idis = (d0Full.norm() != 0),
       reduce_ivel = (v0Full.norm() != 0);
  if(domain->solInfo().useMassNormalizedBasis || domain->solInfo().newmarkBeta == 0) {
    SubDOp *M = NULL;
    if(reduce_idis || reduce_ivel) {
      SparseMatrix **subM = new SparseMatrix * [decDomain->getNumSub()];
      execParal1R(decDomain->getNumSub(), this, &DistrElementSamplingDriver::subMakeMass, subM);
      M = new SubDOp(decDomain->getNumSub(), subM);
    }
    if(reduce_idis) {
      M->mult(d0Full, tmp);
      assembler->assemble(tmp);
      reduce(podBasis, tmp, d0Red);
    }
    if(reduce_ivel) {
      M->mult(v0Full, tmp);
      assembler->assemble(tmp);
      reduce(podBasis, tmp, v0Red);
    }
    if(M) delete M;
  }
  else {
    if(reduce_idis) reduce(podBasis, d0Full, d0Red);
    if(reduce_ivel) reduce(podBasis, v0Full, v0Red);
  }

  if(myID == 0) {
    // Weights output file generation
    outputFullWeights(gweights, gelemIds);

    // Mesh output file generation
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

    // Read in the truncated basis into a (non-distributed) VecBasis XXX
    domain->preProcessing();
    buildDomainCdsa();
    std::string fileName2 = BasisFileId(fileInfo, BasisId::STATE, BasisId::POD);
    if(domain->solInfo().newmarkBeta == 0 || domain->solInfo().useMassNormalizedBasis) fileName2.append(".normalized");
    const VecNodeDof6Conversion vecDofConversion(*domain->getCDSA());
    BasisInputStream in(fileName2, vecDofConversion);
    VecBasis podBasis;
    const int podSizeMax = domain->solInfo().maxSizePodRom;
    if(podSizeMax != 0) {
      readVectors(in, podBasis, podSizeMax);
    } else {
      readVectors(in, podBasis);
    }

    // Output the reduced mesh
    Elemset &inputElemSet = *(geoSource->getElemSet());
    std::auto_ptr<Connectivity> elemToNode(new Connectivity(&inputElemSet));
    const MeshRenumbering meshRenumbering(reducedelemIds.begin(), reducedelemIds.end(), *elemToNode, verboseFlag);
    const MeshDesc reducedMesh(domain, geoSource, meshRenumbering, weightsMap); 
    outputMeshFile(fileInfo, reducedMesh, podBasis.vectorCount());

    // Output the reduced forces
    std::ofstream meshOut(getMeshFilename(fileInfo).c_str(), std::ios_base::app);
    if(domain->solInfo().reduceFollower) meshOut << "EXTFOL\n";
    if(domain->gravityFlag()) {
      meshOut << "*\nFORCES -1\nMODAL\n"; // note: gravity forces are put in loadset -1 so that MFTT (if present) will not be applied
      meshOut.precision(std::numeric_limits<double>::digits10+1);
      for(int i = 0; i < podBasis.vectorCount(); ++i)
        meshOut << i+1 << " " << gravForceRed[i] << std::endl;
    }
    if(reduce_f) {
      meshOut << "*\nFORCES\nMODAL\n";
      meshOut.precision(std::numeric_limits<double>::digits10+1);
      for(int i = 0; i < podBasis.vectorCount(); ++i)
        meshOut << i+1 << " " << constForceRed[i] << std::endl;
    }

    // output the reduced initial conditions
    if(reduce_idis) {
      meshOut << "*\nIDISPLACEMENTS\nMODAL\n";
      meshOut.precision(std::numeric_limits<double>::digits10+1);
      for(int i=0; i<podBasis.vectorCount(); ++i)
        meshOut << i+1 << " " << d0Red[i] << std::endl;
    }
    if(reduce_ivel) {
      meshOut << "*\nIVELOCITIES\nMODAL\n";
      meshOut.precision(std::numeric_limits<double>::digits10+1);
      for(int i=0; i<podBasis.vectorCount(); ++i)
        meshOut << i+1 << " " << v0Red[i] << std::endl;
    }

    meshOut.close();

#ifdef USE_EIGEN3
    // Build and output compressed basis
    DofSetArray reduced_dsa(reducedMesh.nodes().size(), const_cast<Elemset&>(reducedMesh.elements()));
    int num_bc = reducedMesh.dirichletBConds().size();
    BCond *bc = (num_bc > 0) ? const_cast<BCond*>(&reducedMesh.dirichletBConds()[0]) : NULL;
    ConstrainedDSA reduced_cdsa(reduced_dsa, num_bc, bc);
    podBasis.makeSparseBasis(meshRenumbering.reducedNodeIds(), domain->getCDSA(), &reduced_cdsa);
    {
      std::string filename = BasisFileId(fileInfo, BasisId::STATE, BasisId::POD);
      filename.append(".reduced");
      if(domain->solInfo().newmarkBeta == 0 || domain->solInfo().useMassNormalizedBasis) filename.append(".normalized");
      filePrint(stderr," ... Writing compressed basis to file %s ...\n", filename.c_str());
      VecNodeDof6Conversion converter(reduced_cdsa);
      BasisOutputStream output(filename, converter, false);

      for (int iVec = 0; iVec < podBasis.vectorCount(); ++iVec) {
        output << podBasis.compressedBasis().col(iVec);
      }
    }
#endif
  }
}

void
DistrElementSamplingDriver::buildDomainCdsa()
{
  const int numdof = domain->numdof();
  SimpleBuffer<int> bc(numdof);
  SimpleBuffer<double> bcx(numdof);

  domain->make_bc(bc.array(), bcx.array());
  domain->make_constrainedDSA(bc.array());
}

void
DistrElementSamplingDriver::subMakeMass(int i, SparseMatrix **subM)
{
  AllOps<double> allOps;
  allOps.M = subM[i] = decDomain->getSubDomain(i)->constructDBSparseMatrix<double>();
  decDomain->getSubDomain(i)->makeSparseOps<double>(allOps, 0, 0, 0);
}

} /* end namespace Rom */

extern Communicator *structCom;

Rom::DriverInterface *distrElementSamplingDriverNew(Domain *domain)
{
  return new Rom::DistrElementSamplingDriver(domain, structCom);
}
