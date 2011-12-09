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
  void currentTimeIs(double t);
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

  double currentTime_;
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
DistrExplicitSnapshotNonLinDynamic::currentTimeIs(double t) {
  snapshotHandler_->currentTimeIs(t);
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
  skipCounter_(0),
  currentTime_(0.0)
{}

void
DistrExplicitSnapshotNonLinDynamic::SnapshotHandler::currentTimeIs(double time) {
  currentTime_ = time;
}

void
DistrExplicitSnapshotNonLinDynamic::SnapshotHandler::snapshotAdd(const DistrVector &state) {
  ++skipCounter_;
  if (skipCounter_ >= parent_->domain->solInfo().skipPodRom) {
    converter_.paddedNodeDof6(state, buffer_);
    outputFile_.stateAdd(buffer_, currentTime_);
    skipCounter_ = 0;
  }
}

} // end namespace Rom
