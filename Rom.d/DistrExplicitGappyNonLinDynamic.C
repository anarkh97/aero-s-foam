#include "DistrExplicitGappyNonLinDynamic.h"

#include "DistrExplicitGappyProjectionSolver.h"

#include "DistrVecBasis.h"
#include "DistrVecBasisOps.h"

#include "FileNameInfo.h"
#include "DistrBasisFile.h"

#include "DistrMasterMapping.h"
#include "DistrVecNodeDof6Conversion.h"
#include "DistrNodeDof6Buffer.h"
#include "DistrDomainUtils.h"
#include "PtrPtrIterAdapter.h"

#include <Driver.d/DecDomain.h>
#include <Math.d/Vector.h>
#include <Feti.d/DistrVector.h>
#include <Utils.d/DistHelper.h>

#include <algorithm>
#include <stdexcept>
#include <memory>
#include <cstddef>

extern Communicator *structCom;

namespace Rom {

DistrExplicitGappyNonLinDynamic::DistrExplicitGappyNonLinDynamic(Domain *domain) :
  MultiDomainDynam(domain)
{}

void
DistrExplicitGappyNonLinDynamic::preProcess() {
  MultiDomainDynam::preProcess();
 
  FileNameInfo fileInfo; 
  DistrBasisInputFile stateBasisFile(BasisFileId(fileInfo, BasisId::STATE, BasisId::GAPPY_POD));
  DistrBasisInputFile forceBasisFile(BasisFileId(fileInfo, BasisId::FORCE, BasisId::GAPPY_POD));

  const int fileBasisRank = std::min(stateBasisFile.stateCount(), forceBasisFile.stateCount());
  
  const int projectionSubspaceSize = domain->solInfo().maxSizePodRom ?
                                     std::min(domain->solInfo().maxSizePodRom, fileBasisRank) :
                                     fileBasisRank;

  filePrint(stderr, "Projection subspace of dimension = %d\n", projectionSubspaceSize);

  reconstructionBasis_.dimensionIs(projectionSubspaceSize, decDomain->solVecInfo());
  projectionBasis_.dimensionIs(projectionSubspaceSize, decDomain->solVecInfo());

  DistrVecNodeDof6Conversion converter(decDomain->getAllSubDomains(), decDomain->getAllSubDomains() + decDomain->getNumSub());
  
  typedef PtrPtrIterAdapter<SubDomain> SubDomIt;
  DistrMasterMapping masterMapping(SubDomIt(decDomain->getAllSubDomains()),
                                   SubDomIt(decDomain->getAllSubDomains() + decDomain->getNumSub()));
  DistrNodeDof6Buffer snapBuffer(masterMapping.localNodeBegin(), masterMapping.localNodeEnd());

  for (DistrVecBasis::iterator it = reconstructionBasis_.begin(),
                               it_end = reconstructionBasis_.end();
                               it != it_end; ++it) {
    assert(stateBasisFile.validCurrentState());
    stateBasisFile.currentStateBuffer(snapBuffer);
    converter.vector(snapBuffer, *it);
    stateBasisFile.currentStateIndexInc();
  }
  
  for (DistrVecBasis::iterator it = projectionBasis_.begin(),
                               it_end = projectionBasis_.end();
                               it != it_end; ++it) {
    assert(forceBasisFile.validCurrentState());
    forceBasisFile.currentStateBuffer(snapBuffer);
    converter.vector(snapBuffer, *it);
    forceBasisFile.currentStateIndexInc();
  }
}

MDDynamMat *
DistrExplicitGappyNonLinDynamic::buildOps(double mCoef, double cCoef, double kCoef) {
  MDDynamMat *result = MultiDomainDynam::buildOps(mCoef, cCoef, kCoef);
  std::auto_ptr<DistrExplicitGappyProjectionSolver> solver(
      new DistrExplicitGappyProjectionSolver(projectionBasis_, reconstructionBasis_));

  delete result->dynMat;
  result->dynMat = solver.release();

  return result;
}

} // end namespace Rom
