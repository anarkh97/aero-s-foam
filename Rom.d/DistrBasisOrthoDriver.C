#include "DistrBasisOrthoDriver.h"

#include "DistrDomainUtils.h"
#include "DistrMasterMapping.h"
#include "DistrNodeDof6Buffer.h"
#include "DistrVecNodeDof6Conversion.h"
#include "PtrPtrIterAdapter.h"

#include "DistrSvdOrthogonalization.h"

#include "FileNameInfo.h"
#include "DistrBasisFile.h"

#include <Utils.d/DistHelper.h>

#include <algorithm>
#include <memory>
#include <cassert>

namespace Rom {

DistrBasisOrthoDriver::DistrBasisOrthoDriver(Domain *domain, Communicator *comm) :
  domain_(domain),
  comm_(comm)
{}

void
DistrBasisOrthoDriver::solve() {
  std::auto_ptr<DecDomain> decDomain(createDecDomain<double>(domain_));
  decDomain->preProcess();
  
  typedef PtrPtrIterAdapter<SubDomain> SubDomIt;
  DistrMasterMapping masterMapping(SubDomIt(decDomain->getAllSubDomains()),
                                   SubDomIt(decDomain->getAllSubDomains() + decDomain->getNumSub()));

  DistrVecNodeDof6Conversion converter(decDomain->getAllSubDomains(),
                                       decDomain->getAllSubDomains() + decDomain->getNumSub());
  
  FileNameInfo fileInfo;
  const BasisId::Type workload = BasisId::STATE; /*domain_->solInfo().gaussNewtonPodRom ? BasisId::RESIDUAL :
                                 (domain_->solInfo().galerkinPodRom ? BasisId::FORCE : BasisId::STATE);*/

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
}

} /* end namespace Rom */

extern Communicator *structCom;

Rom::DriverInterface *distrBasisOrthoDriverNew(Domain *domain) {
  return new Rom::DistrBasisOrthoDriver(domain, structCom);
}
