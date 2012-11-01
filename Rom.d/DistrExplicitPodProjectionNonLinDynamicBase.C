#include "DistrExplicitPodProjectionNonLinDynamicBase.h"

#include "DistrGalerkinProjectionSolver.h"

#include "DistrVecBasis.h"

#include "FileNameInfo.h"
#include "DistrBasisFile.h"

#include "DistrMasterMapping.h"
#include "DistrVecNodeDof6Conversion.h"
#include "DistrNodeDof6Buffer.h"

#include "PtrPtrIterAdapter.h"

#include <Driver.d/DecDomain.h>
#include <Feti.d/DistrVector.h>
#include <Utils.d/DistHelper.h>

#include <algorithm>
#include <memory>

namespace Rom {

DistrExplicitPodProjectionNonLinDynamicBase::DistrExplicitPodProjectionNonLinDynamicBase(Domain *domain) :
  MultiDomainDynam(domain)
{}

void
DistrExplicitPodProjectionNonLinDynamicBase::preProcess() {
  MultiDomainDynam::preProcess();

  FileNameInfo fileInfo; 
  DistrBasisInputFile podBasisFile(BasisFileId(fileInfo, BasisId::STATE, BasisId::POD));

  const int projectionSubspaceSize = domain->solInfo().maxSizePodRom ?
                                     std::min(domain->solInfo().maxSizePodRom, podBasisFile.stateCount()) :
                                     podBasisFile.stateCount();

  filePrint(stderr, " ... Projection subspace of dimension = %d ...\n", projectionSubspaceSize);
  projectionBasis_.dimensionIs(projectionSubspaceSize, decDomain->masterSolVecInfo());

  DistrVecNodeDof6Conversion converter(decDomain->getAllSubDomains(), decDomain->getAllSubDomains() + decDomain->getNumSub());
  
  typedef PtrPtrIterAdapter<SubDomain> SubDomIt;
  DistrMasterMapping masterMapping(SubDomIt(decDomain->getAllSubDomains()),
                                   SubDomIt(decDomain->getAllSubDomains() + decDomain->getNumSub()));
  DistrNodeDof6Buffer buffer(masterMapping.localNodeBegin(), masterMapping.localNodeEnd());

  for (DistrVecBasis::iterator it = projectionBasis_.begin(),
                               it_end = projectionBasis_.end();
                               it != it_end; ++it) {
    assert(podBasisFile.validCurrentState());

    podBasisFile.currentStateBuffer(buffer);
    converter.vector(buffer, *it);
    
    podBasisFile.currentStateIndexInc();
  }
}

MDDynamMat *
DistrExplicitPodProjectionNonLinDynamicBase::buildOps(double mCoef, double cCoef, double kCoef) {
  MDDynamMat *result = MultiDomainDynam::buildOps(mCoef, cCoef, kCoef);
  assert(result->M);

  const GenSubDOp<double> &fullMass = *(result->M);
  std::auto_ptr<DistrGalerkinProjectionSolver> solver(new DistrGalerkinProjectionSolver(fullMass));

  solver->projectionBasisIs(projectionBasis_);

  delete result->dynMat;
  result->dynMat = solver.release();

  return result;
}

} // end namespace Rom
