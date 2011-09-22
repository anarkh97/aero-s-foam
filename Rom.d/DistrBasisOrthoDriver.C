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
  DistrBasisInputFile inputFile(BasisFileId(fileInfo, BasisId::STATE, BasisId::SNAPSHOTS));

  DistrSvdOrthogonalization solver(comm_, comm_->numCPUs(), 1);
  
  const int blockSize = 1; // TODO more efficient
  {
    solver.blockSizeIs(blockSize);
  }
 
  const int localLength = decDomain->solVecInfo().totLen();
  {
    const int maxLocalLength = comm_->globalMax(localLength);
    const int maxCpuLoad = ((maxLocalLength / blockSize) + (maxLocalLength % blockSize)) * blockSize;
    assert(maxCpuLoad >= localLength);
    const int globalProbSize = maxCpuLoad * solver.rowCpus();
    solver.problemSizeIs(globalProbSize, inputFile.stateCount());
    assert(solver.localRows() == maxCpuLoad);
  }

  {
    DistrNodeDof6Buffer inputBuffer(masterMapping.localNodeBegin(), masterMapping.localNodeEnd());
    while (inputFile.validCurrentState()) {
      inputFile.currentStateBuffer(inputBuffer);
      double *vecBuffer = solver.matrixColBuffer(inputFile.currentStateIndex());
      GenStackDistVector<double> vec(decDomain->solVecInfo(), vecBuffer);
      converter.paddedMasterVector(inputBuffer, vec);
      std::fill(vecBuffer + localLength, vecBuffer + solver.localRows(), 0.0);
      inputFile.currentStateIndexInc();
    }
  }
 
  solver.solve();

  const int podVectorCount = solver.singularValueCount(); // TODO truncate
  {
    DistrBasisOutputFile outputFile(BasisFileId(fileInfo, BasisId::STATE, BasisId::POD), inputFile.nodeCount(), comm_);
    DistrNodeDof6Buffer outputBuffer(masterMapping.masterNodeBegin(), masterMapping.masterNodeEnd());

    for (int iVec = 0; iVec < podVectorCount; ++iVec) {
      double * const vecBuffer = const_cast<double *>(solver.basisColBuffer(iVec));
      const GenStackDistVector<double> vec(decDomain->solVecInfo(), vecBuffer);
      converter.nodeDof6(vec, outputBuffer);
      outputFile.stateAdd(outputBuffer, solver.singularValue(iVec));
    }
  }
}

} /* end namespace Rom */

extern Communicator *structCom;

Rom::DriverInterface *distrBasisOrthoDriverNew(Domain *domain) {
  return new Rom::DistrBasisOrthoDriver(domain, structCom);
}
