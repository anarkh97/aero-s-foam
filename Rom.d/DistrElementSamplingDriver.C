#include "DistrElementSamplingDriver.h"
#include "SubElementSamplingDriver.h"

#include "DistrDomainUtils.h"

#include "DistrVecBasis.h"

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

namespace Rom {

DistrElementSamplingDriver::DistrElementSamplingDriver(Domain *domain, Communicator *comm) :
  domain_(domain),
  comm_(comm)
{}

void
DistrElementSamplingDriver::solve() {
  std::auto_ptr<DecDomain> decDomain(createDecDomain<double>(domain_));
  decDomain->preProcess();
  
  std::cerr << "here in DistrElementSamplingDriver::solve() \n";

  for(int i=0; i<decDomain->getNumSub(); ++i) {

    std::auto_ptr<Rom::DriverInterface> subDriver;
    subDriver.reset(subElementSamplingDriverNew(decDomain->getAllSubDomains()[i]));
    //subDriver->solve();
  }


/*
  typedef PtrPtrIterAdapter<SubDomain> SubDomIt;
  DistrMasterMapping masterMapping(SubDomIt(decDomain->getAllSubDomains()),
                                   SubDomIt(decDomain->getAllSubDomains() + decDomain->getNumSub()));

  DistrVecNodeDof6Conversion converter(decDomain->getAllSubDomains(),
                                       decDomain->getAllSubDomains() + decDomain->getNumSub());
  
  FileNameInfo fileInfo;
  const BasisId::Type workload = BasisId::STATE;

  DistrBasisInputFile inputFile(BasisFileId(fileInfo, workload, BasisId::SNAPSHOTS));

  DistrSvdOrthogonalization solver(comm_, comm_->numCPUs(), 1);
  
  const int blockSize = 64; // TODO More flexible
  {
    solver.blockSizeIs(blockSize);
  }

  const int skipFactor = domain->solInfo().skipPodRom;
  const int basisStateCount = 1 + (inputFile.stateCount() - 1) / skipFactor;
  const int localLength = decDomain->solVecInfo().totLen();
  {
    const int maxLocalLength = comm_->globalMax(localLength);
    const int maxCpuLoad = ((maxLocalLength / blockSize) + (maxLocalLength % blockSize)) * blockSize;
    assert(maxCpuLoad >= localLength);
    const int globalProbSize = maxCpuLoad * solver.rowCpus();
    solver.problemSizeIs(globalProbSize, basisStateCount);
    assert(solver.localRows() == maxCpuLoad);
  }

  {
    DistrNodeDof6Buffer inputBuffer(masterMapping.localNodeBegin(), masterMapping.localNodeEnd());
    
    int count = 0;
    int skipCounter = skipFactor;
    while (count < basisStateCount) {
      assert(inputFile.validCurrentState());
      inputFile.currentStateBuffer(inputBuffer);

      if (skipCounter >= skipFactor) {
        double *vecBuffer = solver.matrixColBuffer(count);
        GenStackDistVector<double> vec(decDomain->solVecInfo(), vecBuffer);
        converter.paddedMasterVector(inputBuffer, vec);
        std::fill(vecBuffer + localLength, vecBuffer + solver.localRows(), 0.0);
        
        skipCounter = 1;
        ++count;
      } else {
        ++skipCounter;
      }

      inputFile.currentStateIndexInc();
    }
  }
*/

  // read basis
  DistrVecBasis podBasis;
  FileNameInfo fileInfo;
  DistrBasisInputFile podBasisFile(BasisFileId(fileInfo, BasisId::STATE, BasisId::POD));

  const int projectionSubspaceSize = domain->solInfo().maxSizePodRom ?
                                     std::min(domain->solInfo().maxSizePodRom, podBasisFile.stateCount()) :
                                     podBasisFile.stateCount();

  filePrint(stderr, " ... Projection subspace of dimension = %d ...\n", projectionSubspaceSize);
  podBasis.dimensionIs(projectionSubspaceSize, decDomain->masterSolVecInfo());

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
  std::cerr << "finished reading basis\n";

  // read snapshots (TODO)
  DistrVecBasis snapshots;
  std::vector<double> timeStamps;
  {
    DistrBasisInputFile in(BasisFileId(fileInfo, BasisId::STATE, BasisId::SNAPSHOTS));
    const int skipFactor = 1; //domain->solInfo().skipPodRom;
    const int skipOffSet = 0; //domain->solInfo().skipOffSet;
    const int basisStateCount = 1 + (in.stateCount() - 1) / skipFactor;
    filePrint(stderr, " ... basisStateCount = %d ...\n", basisStateCount);

    snapshots.dimensionIs(basisStateCount, decDomain->masterSolVecInfo());
    timeStamps.reserve(basisStateCount);

    for (DistrVecBasis::iterator it = snapshots.begin(),
                                 it_end = snapshots.end();
                                 it != it_end; ++it) {
      assert(in.validCurrentState());

      in.currentStateBuffer(buffer);
      converter.vector(buffer, *it);

      in.currentStateIndexInc();
    }

  }
  std::cerr << "finished reading snapshots\n";



  exit(-1);

  for(int i=0; i<decDomain->getNumSub(); ++i) {

    std::auto_ptr<Rom::DriverInterface> subDriver;
    subDriver.reset(subElementSamplingDriverNew(decDomain->getAllSubDomains()[i]));
    subDriver->solve();
  }
 
/*
  solver.solve();

  const int podVectorCount = domain_->solInfo().maxSizePodRom ?
                             std::min(domain_->solInfo().maxSizePodRom, solver.singularValueCount()) :
                             solver.singularValueCount();
  {
    DistrNodeDof6Buffer outputBuffer(masterMapping.masterNodeBegin(), masterMapping.masterNodeEnd());
    DistrBasisOutputFile outputFile(BasisFileId(fileInfo, workload, BasisId::POD),
                                    inputFile.nodeCount(), outputBuffer.globalNodeIndexBegin(), outputBuffer.globalNodeIndexEnd(),
                                    comm_);

    for (int iVec = 0; iVec < podVectorCount; ++iVec) {
      double * const vecBuffer = const_cast<double *>(solver.basisColBuffer(iVec));
      const GenStackDistVector<double> vec(decDomain->solVecInfo(), vecBuffer);
      converter.paddedNodeDof6(vec, outputBuffer);
      outputFile.stateAdd(outputBuffer, solver.singularValue(iVec));
    }
  }
*/
}

} /* end namespace Rom */

extern Communicator *structCom;

Rom::DriverInterface *distrElementSamplingDriverNew(Domain *domain) {
  return new Rom::DistrElementSamplingDriver(domain, structCom);
}
