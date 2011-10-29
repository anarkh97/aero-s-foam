#include "DistrExplicitPodProjectionNonLinDynamic.h"

#include "DistrGalerkinProjectionSolver.h"

#include "DistrVecBasis.h"
#include "DistrVecBasisOps.h"

#include "FileNameInfo.h"
#include "DistrBasisFile.h"

#include "DistrMasterMapping.h"
#include "DistrVecNodeDof6Conversion.h"
#include "DistrNodeDof6Buffer.h"
#include "PtrPtrIterAdapter.h"

#include <Driver.d/DecDomain.h>
#include <Math.d/Vector.h>
#include <Feti.d/DistrVector.h>
#include <Utils.d/DistHelper.h>

#include <Driver.d/GeoSource.h>

#include <algorithm>
#include <stdexcept>
#include <memory>
#include <cstddef>

extern Communicator *structCom;
extern GeoSource *geoSource;

namespace Rom {

class DistrExplicitPodProjectionNonLinDynamic::SnapshotHandler {
public:
  void forceSnapshotAdd(const DistrVector &s);
  explicit SnapshotHandler(DistrExplicitPodProjectionNonLinDynamic *parent);

private:
  typedef PtrPtrIterAdapter<SubDomain> SubDomIt;
  
  DistrExplicitPodProjectionNonLinDynamic *parent_;

  DistrVecNodeDof6Conversion converter_;
  DistrMasterMapping masterMapping_;
  DistrVector assembledSnapshot_;
  DistrNodeDof6Buffer buffer_;
  DistrBasisOutputFile outputFile_;
  
  int skipCounter_;
};

DistrExplicitPodProjectionNonLinDynamic::DistrExplicitPodProjectionNonLinDynamic(Domain *domain) :
  MultiDomainDynam(domain),
  snapshotHandler_(NULL)
{}

DistrExplicitPodProjectionNonLinDynamic::~DistrExplicitPodProjectionNonLinDynamic() {
  delete snapshotHandler_;
}

void
DistrExplicitPodProjectionNonLinDynamic::preProcess() {
  MultiDomainDynam::preProcess();
 
  FileNameInfo fileInfo; 
  DistrBasisInputFile podBasisFile(BasisFileId(fileInfo, BasisId::STATE, BasisId::POD));

  const int projectionSubspaceSize = domain->solInfo().maxSizePodRom ?
                                     std::min(domain->solInfo().maxSizePodRom, podBasisFile.stateCount()) :
                                     podBasisFile.stateCount();

  filePrint(stderr, "Projection subspace of dimension = %d\n", projectionSubspaceSize);
  projectionBasis_.dimensionIs(projectionSubspaceSize, decDomain->masterSolVecInfo());

  DistrVecNodeDof6Conversion converter(decDomain->getAllSubDomains(), decDomain->getAllSubDomains() + decDomain->getNumSub());
  
  typedef PtrPtrIterAdapter<SubDomain> SubDomIt;
  DistrMasterMapping masterMapping(SubDomIt(decDomain->getAllSubDomains()),
                                   SubDomIt(decDomain->getAllSubDomains() + decDomain->getNumSub()));
  DistrNodeDof6Buffer snapBuffer(masterMapping.localNodeBegin(), masterMapping.localNodeEnd());

  for (DistrVecBasis::iterator it = projectionBasis_.begin(),
                               it_end = projectionBasis_.end();
                               it != it_end; ++it) {
    assert(podBasisFile.validCurrentState());

    podBasisFile.currentStateBuffer(snapBuffer);
    converter.vector(snapBuffer, *it);
    
    podBasisFile.currentStateIndexInc();
  }

  snapshotHandler_ = new SnapshotHandler(this);
}

MDDynamMat *
DistrExplicitPodProjectionNonLinDynamic::buildOps(double mCoef, double cCoef, double kCoef) {
  MDDynamMat *result = MultiDomainDynam::buildOps(mCoef, cCoef, kCoef);
  assert(result->M);

  const GenSubDOp<double> &fullMass = *(result->M);
  std::auto_ptr<DistrGalerkinProjectionSolver> solver(new DistrGalerkinProjectionSolver(fullMass));
  solver->projectionBasisIs(projectionBasis_);

  delete result->dynMat;
  result->dynMat = solver.release();

  return result;
}

void
DistrExplicitPodProjectionNonLinDynamic::forceSnapshotAdd(const DistrVector &f) {
  snapshotHandler_->forceSnapshotAdd(f);
}

DistrExplicitPodProjectionNonLinDynamic::SnapshotHandler::SnapshotHandler(DistrExplicitPodProjectionNonLinDynamic *parent) :
  parent_(parent),
  converter_(parent->decDomain->getAllSubDomains(), parent->decDomain->getAllSubDomains() + parent->decDomain->getNumSub()),
  masterMapping_(SubDomIt(parent->decDomain->getAllSubDomains()), SubDomIt(parent->decDomain->getAllSubDomains() + parent->decDomain->getNumSub())),
  buffer_(masterMapping_.masterNodeBegin(), masterMapping_.masterNodeEnd()),
  assembledSnapshot_(parent->decDomain->solVecInfo()),
  outputFile_(BasisFileId(FileNameInfo(), BasisId::FORCE, BasisId::SNAPSHOTS), geoSource->getNumGlobNodes(),
              buffer_.globalNodeIndexBegin(), buffer_.globalNodeIndexEnd(), structCom),
  skipCounter_(0)
{}

void
DistrExplicitPodProjectionNonLinDynamic::SnapshotHandler::forceSnapshotAdd(const DistrVector &f) {
  ++skipCounter_;
  if (skipCounter_ >= parent_->domain->solInfo().skipPodRom) {
    assembledSnapshot_ = f;
    parent_->decDomain->getSolVecAssembler()->assemble(assembledSnapshot_);
    converter_.paddedNodeDof6(assembledSnapshot_, buffer_);
    outputFile_.stateAdd(buffer_);
    skipCounter_ = 0;
  }
}

} // end namespace Rom
