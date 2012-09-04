#include "DistrSnapshotNonLinDynamic.h"

#include "DistrBasisFile.h"
#include "FileNameInfo.h"
#include "DistrNodeDof6Buffer.h"
#include "DistrMasterMapping.h"
#include "PtrPtrIterAdapter.h"

#include <Driver.d/Domain.h>
#include <Driver.d/DecDomain.h>
#include <Driver.d/GeoSource.h>

#include <Utils.d/dofset.h>
#include <Corotational.d/utilities.h>

#include <cstddef>
#include <memory>

extern Communicator *structCom;
extern GeoSource *geoSource;

namespace Rom {

// Dummy class holding the implementation of DistrSnapshotNonLinDynamic
struct DistrSnapshotNonLinDynamicDetail : private DistrSnapshotNonLinDynamic {
  class RawImpl : public Impl {
  public:
    // Overriden
    virtual void lastMidTimeIs(double t);
    virtual void lastDeltaIs(double dt);
    virtual void stateSnapshotAdd(const DistrGeomState &);
    virtual void postProcess();

    int nodeCount() const { return geoSource->getNumGlobNodes(); }

    explicit RawImpl(DecDomain *);
  
  private:
    typedef PtrPtrIterAdapter<SubDomain> SubDomIt;
    
    DecDomain * decDomain_;
    DistrMasterMapping masterMapping_;

    DistrNodeDof6Buffer snapBuffer_;
    
    FileNameInfo fileInfo_;
    DistrBasisOutputFile stateSnapFile_;
    
    double timeStamp_;
    int stateSkip_;
  };

private:
  // Dummy constructor to avoid compilation failures
  DistrSnapshotNonLinDynamicDetail() :
    DistrSnapshotNonLinDynamic(NULL)
  {}
};

DistrSnapshotNonLinDynamicDetail::RawImpl::RawImpl(DecDomain *decDomain) :
  decDomain_(decDomain),
  masterMapping_(SubDomIt(decDomain->getAllSubDomains()), SubDomIt(decDomain->getAllSubDomains() + decDomain->getNumSub())),
  snapBuffer_(masterMapping_.masterNodeBegin(), masterMapping_.masterNodeEnd()),
  fileInfo_(),
  stateSnapFile_(BasisFileId(fileInfo_, BasisId::STATE, BasisId::SNAPSHOTS), nodeCount(),
                 snapBuffer_.globalNodeIndexBegin(), snapBuffer_.globalNodeIndexEnd(), structCom),
  timeStamp_(decDomain->getDomain()->solInfo().initialTime),
  stateSkip_(0)
{}

void
DistrSnapshotNonLinDynamicDetail::RawImpl::postProcess() {
  // Nothing to do
}

void
DistrSnapshotNonLinDynamicDetail::RawImpl::lastMidTimeIs(double t) {
  timeStamp_ = t;
}

void
DistrSnapshotNonLinDynamicDetail::RawImpl::lastDeltaIs(double dt) {
  timeStamp_ += dt;
}

void
DistrSnapshotNonLinDynamicDetail::RawImpl::stateSnapshotAdd(const DistrGeomState &snap) {
  ++stateSkip_;
  if(stateSkip_ >= decDomain_->getDomain()->solInfo().skipState) {
  const int subDomCount = snap.getNumSub();
  DistrMasterMapping::SubMasterMappingIt mappingIt = masterMapping_.begin();
  for (int iSub = 0; iSub < subDomCount; ++iSub) {
    const GeomState &subSnap = *snap[iSub];
    const CoordSet &refCoords = decDomain_->getSubDomain(iSub)->getNodes();
    const MasterMapping &mapping = *mappingIt++;

    typedef MasterMapping::IndexPairIterator IndexPairIt;
    const IndexPairIt nodeItEnd = mapping.end();
    for (IndexPairIt nodeIt = mapping.begin(); nodeIt != nodeItEnd; ++nodeIt) {
      // Indexing
      const int iLocalNode = nodeIt->local;
      const int iGlobalNode = nodeIt->global;

      // Translational dofs
      snapBuffer_[iGlobalNode][0] = subSnap[iLocalNode].x - refCoords[iLocalNode]->x;
      snapBuffer_[iGlobalNode][1] = subSnap[iLocalNode].y - refCoords[iLocalNode]->y;
      snapBuffer_[iGlobalNode][2] = subSnap[iLocalNode].z - refCoords[iLocalNode]->z;

      // Rotational dofs
      mat_to_vec(const_cast<double (*)[3]>(subSnap[iLocalNode].R), &snapBuffer_[iGlobalNode][3]);
    }
  }

  stateSnapFile_.stateAdd(snapBuffer_, timeStamp_);
  stateSkip_ = 0;
 }
}


DistrSnapshotNonLinDynamic::DistrSnapshotNonLinDynamic(Domain *domain) :
  MDNLDynamic(domain),
  impl_(NULL)
{}

void
DistrSnapshotNonLinDynamic::preProcess() {
  MDNLDynamic::preProcess();
  impl_.reset(new DistrSnapshotNonLinDynamicDetail::RawImpl(this->getDecDomain()));
}

} /* end namespace Rom */
