#include "DistrSnapshotClusteringDriver.h"

#include <Utils.d/DistHelper.h>

#if defined (USE_SCALAPACK) && defined (USE_EIGEN3)
#include "DistrDomainUtils.h"
#include "DistrMasterMapping.h"
#include "DistrNodeDof6Buffer.h"
#include "DistrVecNodeDof6Conversion.h"
#include "PtrPtrIterAdapter.h"
#include "DistrVecBasis.h"
#include "DistrVecBasisOps.h"

#include "DistrSnapshotClusteringSolver.h"

#include "FileNameInfo.h"
#include "DistrBasisFile.h"

#include <algorithm>
#include <memory>
#include <cassert>
#include <sstream>

namespace Rom {

DistrSnapshotClusteringDriver::DistrSnapshotClusteringDriver(Domain *domain, Communicator *comm) :
  MultiDomainDynam(domain),
  domain_(domain),
  comm_(comm)
{}

//Non-member functions
//===============
void readIntoSolver(DistrSnapshotClusteringSolver &solver, DistrNodeDof6Buffer &inputBuffer, int& nodeCount,
                    DistrVecNodeDof6Conversion &converter, const DistrInfo &vectorSize, BasisId::Type workload,
                    BasisId::Level type, int numEntries, int& solverCol, int skipFactor=1)
{
  FileNameInfo fileInfo;
  for(int i = 0; i < numEntries; i++) {
    std::string fileName = BasisFileId(fileInfo,workload,type,i);
    DistrBasisInputFileTemplate<6> inputFile(fileName);
    filePrint(stderr, " ... Reading in Snapshot file: %s ...\n", fileName.c_str());
    nodeCount = inputFile.nodeCount();
    int basisStateCount = 1+(inputFile.stateCount()-1)/skipFactor;
    {
      int count = 0;
      int skipCounter = skipFactor;
      while (count < basisStateCount) {
        assert(inputFile.validCurrentState());
        inputFile.currentStateBuffer(inputBuffer);

        if (skipCounter >= skipFactor) {
          double *vecBuffer = solver.matrixColBuffer(solverCol);
          GenStackDistVector<double> vec(vectorSize, vecBuffer);

          converter.unpaddedMasterVector(inputBuffer, vec);
          std::fill(vecBuffer + vectorSize.totLen(), vecBuffer + solver.localRows(), 0.0);
          skipCounter = 1;
          ++solverCol;
          ++count;
        } else {
          ++skipCounter;
        }
        inputFile.currentStateIndexInc();
      }
    }
  }
}

void writeOutofSolver(DistrSnapshotClusteringSolver &solver, DistrNodeDof6Buffer &outputBuffer, int nodeCount,
                      DistrVecNodeDof6Conversion &converter, const DistrInfo &vectorSize, BasisId::Type workload,
                      BasisId::Level type, int numClusters, Communicator *comm)
{
  FileNameInfo fileInfo;
  for(int i=0; i<numClusters; ++i) {
    int clusterDim = solver.clusterColCount(i);
    std::string fileName = BasisFileId(fileInfo, workload, type);
    std::ostringstream ss;
    ss << ".cluster" << i+1;
    fileName.append(ss.str());
    DistrBasisOutputFile outputFile(fileName,
                                    nodeCount, outputBuffer.globalNodeIndexBegin(), outputBuffer.globalNodeIndexEnd(),
                                    comm, false);
    filePrint(stderr, " ... Writing %d clustered snapshots to file %s ...\n", clusterDim, fileName.c_str());
    for (int iVec = 0; iVec < clusterDim; ++iVec) {
      double * const vecBuffer = const_cast<double *>(solver.clusterColBuffer(i,iVec));
      const GenStackDistVector<double> vec(vectorSize, vecBuffer);
      converter.unpaddedNodeDof6(vec, outputBuffer);
      outputFile.stateAdd(outputBuffer, 1.0);
    }
    fileName.append(".centroid");
    DistrBasisOutputFile outputFile2(fileName,
                                     nodeCount, outputBuffer.globalNodeIndexBegin(), outputBuffer.globalNodeIndexEnd(),
                                     comm, false);
    filePrint(stderr, " ... Writing centroid of %d clustered snapshots to file %s ...\n", clusterDim, fileName.c_str());
    {
      double * const vecBuffer = const_cast<double *>(solver.clusterCentroidBuffer(i));
      const GenStackDistVector<double> vec(vectorSize, vecBuffer);
      converter.unpaddedNodeDof6(vec, outputBuffer);
      outputFile2.stateAdd(outputBuffer, 1.0);
    }
  }
}

void
DistrSnapshotClusteringDriver::solve() {
  
  MultiDomainDynam::preProcess();

  const int blockSize = domain->solInfo().svdBlockSize;

  typedef PtrPtrIterAdapter<SubDomain> SubDomIt;
  for(int i=0; i<decDomain->getNumSub(); ++i) decDomain->getSubDomain(i)->computeMasterFlag(decDomain->getMpcToSub());
  DistrMasterMapping masterMapping(SubDomIt(decDomain->getAllSubDomains()),
                                   SubDomIt(decDomain->getAllSubDomains() + decDomain->getNumSub()));

  DistrVecNodeDof6Conversion converter(decDomain->getAllSubDomains(),
                                       decDomain->getAllSubDomains() + decDomain->getNumSub());

  FileNameInfo fileInfo;
  const BasisId::Type workload = BasisId::STATE;

  const int skipFactor = domain->solInfo().skipPodRom;
  int stateCount = 0;
  int nodeCount = 0;
  int snapBasisStateCount = 0;
  if(!domain->solInfo().snapfiPodRom.empty()) {
    for(int i = 0; i < domain->solInfo().snapfiPodRom.size(); i++) {
      std::string fileName = BasisFileId(fileInfo, workload, BasisId::SNAPSHOTS, i);
      DistrBasisInputFileTemplate<6> inputFile(fileName);
      snapBasisStateCount += 1+(inputFile.stateCount()-1)/skipFactor;
    }
  }
  else {
    filePrint(stderr, " *** ERROR: no files provided\n");
    exit(-1);
  }

  DistrInfo distrInfo; decDomain->makeNonOverlappingDistrInfo(distrInfo);

  int localLength = distrInfo.totLen();
  const int numClusters = domain->solInfo().clustering;
  const int globalProbSize = domain->getCDSA()->size();
 
  DistrSnapshotClusteringSolver solver(comm_, globalProbSize, snapBasisStateCount, localLength, numClusters, blockSize);
 
  int solverCol = 0;
  DistrNodeDof6Buffer inputBuffer(masterMapping.localNodeBegin(), masterMapping.localNodeEnd());
  readIntoSolver(solver, inputBuffer, nodeCount, converter, distrInfo, BasisId::STATE, BasisId::SNAPSHOTS,
                 domain->solInfo().snapfiPodRom.size(), solverCol, domain->solInfo().skipPodRom); // read in snapshots

  filePrint(stderr, " ... Partitioning snapshots into %d clusters ...\n", numClusters);
  solver.solve();

  DistrNodeDof6Buffer outputBuffer(masterMapping.masterNodeBegin(), masterMapping.masterNodeEnd());
  writeOutofSolver(solver, outputBuffer, nodeCount, converter, distrInfo, BasisId::STATE, BasisId::SNAPSHOTS,
                   numClusters, comm_);

  if(!domain_->solInfo().velocPodRomFile.empty()) {
    solverCol = 0;
    readIntoSolver(solver, inputBuffer, nodeCount, converter, distrInfo, BasisId::VELOCITY, BasisId::SNAPSHOTS,
                   domain->solInfo().velocPodRomFile.size(), solverCol, domain->solInfo().skipPodRom); // read in velocity snapshots
    writeOutofSolver(solver, outputBuffer, nodeCount, converter, distrInfo, BasisId::VELOCITY, BasisId::SNAPSHOTS,
                     numClusters, comm_);
  }

  if(!domain_->solInfo().accelPodRomFile.empty()) {
    solverCol = 0;
    readIntoSolver(solver, inputBuffer, nodeCount, converter, distrInfo, BasisId::ACCELERATION, BasisId::SNAPSHOTS,
                   domain->solInfo().accelPodRomFile.size(), solverCol, domain->solInfo().skipPodRom); // read in acceleration snapshots
    writeOutofSolver(solver, outputBuffer, nodeCount, converter, distrInfo, BasisId::ACCELERATION, BasisId::SNAPSHOTS,
                     numClusters, comm_);
  }
}

} /* end namespace Rom */

extern Communicator *structCom;

Rom::DriverInterface *distrSnapshotClusteringDriverNew(Domain *domain) {
  return new Rom::DistrSnapshotClusteringDriver(domain, structCom);
}

#else

Rom::DriverInterface *distrSnapshotClusteringDriverNew(Domain *domain) {
  filePrint(stderr, " *** ERROR: requested driver requires ScaLAPACK and Eigen libraries\n");
  exit(-1);
  return 0;
}

#endif
