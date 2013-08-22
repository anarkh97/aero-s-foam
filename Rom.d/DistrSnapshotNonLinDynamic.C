#include "DistrSnapshotNonLinDynamic.h"

#include "DistrBasisFile.h"
#include "FileNameInfo.h"
#include "DistrNodeDof6Buffer.h"
#include "DistrMasterMapping.h"
#include "DistrVecNodeDof6Conversion.h"
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
    virtual void velocSnapshotAdd(const DistrVector &);
    virtual void accelSnapshotAdd(const DistrVector &);
    virtual void postProcess();

    int nodeCount() const { return geoSource->getNumGlobNodes(); }

    explicit RawImpl(DecDomain *);
  
  private:
    typedef PtrPtrIterAdapter<SubDomain> SubDomIt;
    
    DecDomain * decDomain_;
    DistrVecNodeDof6Conversion converter_;
    DistrMasterMapping masterMapping_;

    DistrNodeDof6Buffer snapBuffer_;
    
    FileNameInfo fileInfo_;
    DistrBasisOutputFile *stateSnapFile_;
    DistrBasisOutputFile *velocSnapFile_;
    DistrBasisOutputFile *accelSnapFile_;
    
    double timeStamp_;
    int stateSkip_;
    int velocSkip_;
    int accelSkip_;
  };

private:
  // Dummy constructor to avoid compilation failures
  DistrSnapshotNonLinDynamicDetail() :
    DistrSnapshotNonLinDynamic(NULL)
  {}
};

DistrSnapshotNonLinDynamicDetail::RawImpl::RawImpl(DecDomain *decDomain) :
  decDomain_(decDomain),
  converter_(decDomain->getAllSubDomains(), decDomain->getAllSubDomains() + decDomain->getNumSub()),
  masterMapping_(SubDomIt(decDomain->getAllSubDomains()), SubDomIt(decDomain->getAllSubDomains() + decDomain->getNumSub())),
  snapBuffer_(masterMapping_.masterNodeBegin(), masterMapping_.masterNodeEnd()),
  fileInfo_(),
  stateSnapFile_(NULL),
  velocSnapFile_(NULL),
  accelSnapFile_(NULL),
  timeStamp_(decDomain->getDomain()->solInfo().initialTime),
  stateSkip_(0),
  velocSkip_(0),
  accelSkip_(0)
{
  if(decDomain->getDomain()->solInfo().statevectPodRom){
    stateSnapFile_ = new DistrBasisOutputFile(BasisFileId(fileInfo_, BasisId::STATE, BasisId::SNAPSHOTS), nodeCount(),
                 snapBuffer_.globalNodeIndexBegin(), snapBuffer_.globalNodeIndexEnd(), structCom,
                 (geoSource->getCheckFileInfo()->lastRestartFile != 0));}
  if(decDomain->getDomain()->solInfo().velocvectPodRom){
    velocSnapFile_ = new DistrBasisOutputFile(BasisFileId(fileInfo_, BasisId::VELOCITY, BasisId::SNAPSHOTS), nodeCount(),
                 snapBuffer_.globalNodeIndexBegin(), snapBuffer_.globalNodeIndexEnd(), structCom,
                 (geoSource->getCheckFileInfo()->lastRestartFile != 0));}
  if(decDomain->getDomain()->solInfo().accelvectPodRom){
    accelSnapFile_ = new DistrBasisOutputFile(BasisFileId(fileInfo_, BasisId::ACCELERATION, BasisId::SNAPSHOTS), nodeCount(),
                 snapBuffer_.globalNodeIndexBegin(), snapBuffer_.globalNodeIndexEnd(), structCom,
                 (geoSource->getCheckFileInfo()->lastRestartFile != 0));}
}

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
  if(stateSnapFile_ && (stateSkip_ >= decDomain_->getDomain()->solInfo().skipState)) {
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
      // old method: collect the rescaled rotation vector
      //mat_to_vec(const_cast<double (*)[3]>(subSnap[iLocalNode].R), &snapBuffer_[iGlobalNode][3]);

      // new method: collect the unscaled rotation vector which has already be computed and stored in NodeState::theta
      snapBuffer_[iGlobalNode][3] = subSnap[iLocalNode].theta[0];
      snapBuffer_[iGlobalNode][4] = subSnap[iLocalNode].theta[1];
      snapBuffer_[iGlobalNode][5] = subSnap[iLocalNode].theta[2];
    }
  }

  stateSnapFile_->stateAdd(snapBuffer_, timeStamp_);
  stateSkip_ = 0;
 }
}

void
DistrSnapshotNonLinDynamicDetail::RawImpl::velocSnapshotAdd(const DistrVector &veloc) {
  ++velocSkip_;
  if (velocSnapFile_ && (velocSkip_ >= decDomain_->getDomain()->solInfo().skipVeloc)) {
    converter_.paddedNodeDof6(veloc, snapBuffer_);
    velocSnapFile_->stateAdd(snapBuffer_, timeStamp_);
    velocSkip_ = 0;
  }
}

void
DistrSnapshotNonLinDynamicDetail::RawImpl::accelSnapshotAdd(const DistrVector &accel) {
  ++accelSkip_;
  if (accelSnapFile_ && (accelSkip_ >= decDomain_->getDomain()->solInfo().skipAccel)) {
    converter_.paddedNodeDof6(accel, snapBuffer_);
    accelSnapFile_->stateAdd(snapBuffer_, timeStamp_);
    accelSkip_ = 0;
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

void
DistrSnapshotNonLinDynamic::saveVelocSnapshot(DistrGeomState &state, const DistrVector &veloc) {
  if(domain->solInfo().velocvectPodRom) {
    if(state.getHaveRot()) {
      DistrVector v(veloc);
      state.transform(v, 2); // transform convected angular velocity to time derivative of total rotation vector
      impl_->velocSnapshotAdd(v);
    }
    else {
      impl_->velocSnapshotAdd(veloc);
    }
  }

}

void
DistrSnapshotNonLinDynamic::saveAccelSnapshot(DistrGeomState &state, const DistrVector &accel) {
  if(domain->solInfo().accelvectPodRom) {
    if(state.getHaveRot()) {
      DistrVector a(accel);
      state.transform(a, 6); // transform convected angular acceleration to second time derivative of total rotation vector
      impl_->accelSnapshotAdd(a);
    }
    else {
      impl_->accelSnapshotAdd(accel);
    }
  }
}

} /* end namespace Rom */
