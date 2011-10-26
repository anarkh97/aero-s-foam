#include "DistrExplicitSnapshotNonLinDynamic.h"

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

#include <Driver.d/GeoSource.h>

#include <algorithm>
#include <cstddef>

extern Communicator *structCom;
extern GeoSource *geoSource;

namespace Rom {

class DistrExplicitSnapshotNonLinDynamic::SnapshotHandler {
public:
  void snapshotAdd(const DistrVector &s);
  explicit SnapshotHandler(DistrExplicitSnapshotNonLinDynamic *parent);

private:
  typedef PtrPtrIterAdapter<SubDomain> SubDomIt;
  
  DistrExplicitSnapshotNonLinDynamic *parent_;

  DistrVecNodeDof6Conversion converter_;
  DistrMasterMapping masterMapping_;
  DistrNodeDof6Buffer buffer_;
  DistrBasisOutputFile outputFile_;
  
  int skipCounter_;
};

DistrExplicitSnapshotNonLinDynamic::DistrExplicitSnapshotNonLinDynamic(Domain *domain) :
  MultiDomainDynam(domain),
  snapshotHandler_(NULL)
{}

DistrExplicitSnapshotNonLinDynamic::~DistrExplicitSnapshotNonLinDynamic() {
  delete snapshotHandler_;
}

void
DistrExplicitSnapshotNonLinDynamic::preProcess() {
  MultiDomainDynam::preProcess();
  snapshotHandler_ = new SnapshotHandler(this);
}

void
DistrExplicitSnapshotNonLinDynamic::snapshotAdd(const DistrVector &f) {
  snapshotHandler_->snapshotAdd(f);
}

DistrExplicitSnapshotNonLinDynamic::SnapshotHandler::SnapshotHandler(DistrExplicitSnapshotNonLinDynamic *parent) :
  parent_(parent),
  converter_(parent->decDomain->getAllSubDomains(), parent->decDomain->getAllSubDomains() + parent->decDomain->getNumSub()),
  masterMapping_(SubDomIt(parent->decDomain->getAllSubDomains()), SubDomIt(parent->decDomain->getAllSubDomains() + parent->decDomain->getNumSub())),
  buffer_(masterMapping_.masterNodeBegin(), masterMapping_.masterNodeEnd()),
  outputFile_(BasisFileId(FileNameInfo(), BasisId::STATE, BasisId::SNAPSHOTS), geoSource->getNumGlobNodes(),
              buffer_.globalNodeIndexBegin(), buffer_.globalNodeIndexEnd(), structCom),
  skipCounter_(0)
{}

void
DistrExplicitSnapshotNonLinDynamic::SnapshotHandler::snapshotAdd(const DistrVector &acc) {
  ++skipCounter_;
  if (skipCounter_ >= parent_->domain->solInfo().skipPodRom) {
    converter_.nodeDof6(acc, buffer_);
    outputFile_.stateAdd(buffer_);
    skipCounter_ = 0;
  }
}

} // end namespace Rom
