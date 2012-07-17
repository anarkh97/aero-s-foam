#include "SnapshotNonLinDynamic.h"

#include "SvdOrthogonalization.h"
#include "VecNodeDof6Conversion.h"
#include "BasisBinaryFile.h"
#include "NodeDof6Buffer.h"
#include "FileNameInfo.h"
#include "BasisFileStream.h"

#include <Driver.d/Domain.h>

#include <Utils.d/dofset.h>
#include <Corotational.d/utilities.h>

#include <deque>

#include <cstddef>
#include <algorithm>
#include <memory>

namespace Rom {

// Dummy class holding the implementation of SnapshotNonLinDynamic
struct SnapshotNonLinDynamicDetail : private SnapshotNonLinDynamic {
  class sttSnapImpl : public Impl {
  public:
     void lastMidTimeIs(double t);
     void lastDeltaIs(double dt);
     void stateSnapshotAdd(const GeomState &);
     void handleResidualSnapshot(const Vector &res);
     void handleJacobianSnapshot();
     void postProcess();

    int dofSetNodeCount() const { return converter_.dofSetNodeCount(); }
    int vectorSize() const { return converter_.vectorSize(); }

    explicit sttSnapImpl(Domain *, BasisId::Level level = BasisId::SNAPSHOTS);

  protected:
    template <typename VecType>
    void fillSnapBuffer(const VecType &origin);
  
    const NodeDof6Buffer &snapBuffer() const { return snapBuffer_; }
    const VecNodeDof6Conversion &converter() const { return converter_; }

    int maxSizePodRom() const { return domain_->solInfo().maxSizePodRom; }

  private:
    Domain * domain_;
    
    VecNodeDof6Conversion converter_;
    NodeDof6Buffer snapBuffer_;

    double timeStamp_;

  protected:
    FileNameInfo fileInfo_;
    BasisBinaryOutputFile stateSnapFile_;
  };

// Implementation with residual snapshots
class resSnapImpl : public Impl {
public:
  explicit resSnapImpl(SnapshotNonLinDynamic *parent, Domain *domain);

  // Overriden functions
   void lastMidTimeIs(double t);
   void lastDeltaIs(double dt);
   void stateSnapshotAdd(const GeomState &state);
   void handleResidualSnapshot(const Vector &res);
   void handleJacobianSnapshot();
   void postProcess();

private:
  SnapshotNonLinDynamic *parent_;
  VecNodeDof6Conversion vecNodeDof6Conversion_;
  FileNameInfo fileInfo_;
  BasisOutputStream residualSnapFile_;
};

//Implementation with jacobian snapshots
class jacSnapImpl : public Impl {
public:
  explicit jacSnapImpl(SnapshotNonLinDynamic *parent, Domain *domain);

  // Overriden functions
   void lastMidTimeIs(double t);
   void lastDeltaIs(double dt);
   void stateSnapshotAdd(const GeomState &state);
   void handleResidualSnapshot(const Vector &res);
   void handleJacobianSnapshot();
   void postProcess();

private:
  SnapshotNonLinDynamic *parent_;
  VecNodeDof6Conversion vecNodeDof6Conversion_;
  FileNameInfo fileInfo_;
  BasisOutputStream jacobianSnapFile_;
};

private:
  // Dummy constructor to avoid compilation failures
  SnapshotNonLinDynamicDetail(Domain *d) :
    SnapshotNonLinDynamic(d)
  {}
};

SnapshotNonLinDynamicDetail::sttSnapImpl::sttSnapImpl(Domain * domain, BasisId::Level level) :
  domain_(domain),
  converter_(*domain->getCDSA()),
  snapBuffer_(dofSetNodeCount()),
  fileInfo_(),
  stateSnapFile_(BasisFileId(fileInfo_, BasisId::STATE, level), dofSetNodeCount()),
  timeStamp_(domain->solInfo().initialTime)
{}

SnapshotNonLinDynamicDetail::resSnapImpl::resSnapImpl(SnapshotNonLinDynamic *parent, Domain *domain) :
  parent_(parent),
  vecNodeDof6Conversion_(*domain->getCDSA()),
  fileInfo_(),
  residualSnapFile_(BasisFileId(fileInfo_, BasisId::RESIDUAL, BasisId::SNAPSHOTS), vecNodeDof6Conversion_)
{}

SnapshotNonLinDynamicDetail::jacSnapImpl::jacSnapImpl(SnapshotNonLinDynamic *parent, Domain * domain) :
  parent_(parent),
  vecNodeDof6Conversion_(*domain->getCDSA()),
  fileInfo_(),
  jacobianSnapFile_(BasisFileId(fileInfo_, BasisId::JACOBIAN, BasisId::SNAPSHOTS), vecNodeDof6Conversion_)
{}

void
SnapshotNonLinDynamicDetail::sttSnapImpl::postProcess() {
  // Nothing to do
}

void
SnapshotNonLinDynamicDetail::resSnapImpl::postProcess() {
  // Nothing to do
}

void
SnapshotNonLinDynamicDetail::jacSnapImpl::postProcess() {
  // Nothing to do
}

template <typename VecType>
inline
void
SnapshotNonLinDynamicDetail::sttSnapImpl::fillSnapBuffer(const VecType &snap) {
  converter_.paddedNodeDof6(snap, snapBuffer_);
}

void
SnapshotNonLinDynamicDetail::sttSnapImpl::lastMidTimeIs(double t) {
  timeStamp_ = t;
}

void
SnapshotNonLinDynamicDetail::sttSnapImpl::lastDeltaIs(double dt) {
  timeStamp_ += dt;
}

void
SnapshotNonLinDynamicDetail::resSnapImpl::lastMidTimeIs(double t) {
  //empty
}

void
SnapshotNonLinDynamicDetail::resSnapImpl::lastDeltaIs(double dt) {
  //empty
}

void
SnapshotNonLinDynamicDetail::jacSnapImpl::lastMidTimeIs(double t) {
  //empty
}

void
SnapshotNonLinDynamicDetail::jacSnapImpl::lastDeltaIs(double dt) {
  //empty
}

void
SnapshotNonLinDynamicDetail::sttSnapImpl::stateSnapshotAdd(const GeomState &snap) {
  const CoordSet &refCoords = domain_->getNodes();

  for (int iNode = 0, iNodeEnd = dofSetNodeCount(); iNode != iNodeEnd; ++iNode) {
    double *nodeBuffer = snapBuffer_[iNode];

    const Node *refNode = refCoords[iNode];
    if (refNode) {
      const NodeState &snapNode = snap[iNode];

      // Translational dofs
      nodeBuffer[0] = snapNode.x - refNode->x;
      nodeBuffer[1] = snapNode.y - refNode->y;
      nodeBuffer[2] = snapNode.z - refNode->z;

      // Rotational dofs
      mat_to_vec(const_cast<double (*)[3]>(snapNode.R), &nodeBuffer[3]);
    } else {
      // Node does not really exist, corresponds to a gap in node numbering
      std::fill_n(nodeBuffer, 6, 0.0);
    }
  }

  stateSnapFile_.stateAdd(snapBuffer_, timeStamp_);
}

void
SnapshotNonLinDynamicDetail::resSnapImpl::stateSnapshotAdd(const GeomState &snap) {}

void
SnapshotNonLinDynamicDetail::jacSnapImpl::stateSnapshotAdd(const GeomState &snap) {}

void
SnapshotNonLinDynamicDetail::sttSnapImpl::handleResidualSnapshot(const Vector &res) {
  //empty
}

void
SnapshotNonLinDynamicDetail::resSnapImpl::handleResidualSnapshot(const Vector &res) {
  residualSnapFile_ << res;
}

void
SnapshotNonLinDynamicDetail::jacSnapImpl::handleResidualSnapshot(const Vector &res) {
  //empty
}

void
SnapshotNonLinDynamicDetail::sttSnapImpl::handleJacobianSnapshot() {
  //empty
}

void
SnapshotNonLinDynamicDetail::resSnapImpl::handleJacobianSnapshot() {
  //empty
}

void
SnapshotNonLinDynamicDetail::jacSnapImpl::handleJacobianSnapshot() {
  Vector snap(parent_->solVecInfo());
//  expand(getSolver()->lastReducedMatrixAction(), getSolver()->lastReducedSolution(), snap);
  jacobianSnapFile_ << snap;
}

SnapshotNonLinDynamic::SnapshotNonLinDynamic(Domain *domain) :
  NonLinDynamic(domain),
  stateImpl_(NULL),
  resImpl_(NULL),
  jacImpl_(NULL)
{}

void
SnapshotNonLinDynamic::preProcess() {
  NonLinDynamic::preProcess();
  if(domain->solInfo().statevectPodRom)
    stateImpl_.reset(new SnapshotNonLinDynamicDetail::sttSnapImpl(this->domain));
  if(domain->solInfo().residvectPodRom)
    resImpl_.reset(new SnapshotNonLinDynamicDetail::resSnapImpl(this,this->domain));
  if(domain->solInfo().jacobvectPodRom)
    jacImpl_.reset(new SnapshotNonLinDynamicDetail::jacSnapImpl(this,this->domain));
}

void
SnapshotNonLinDynamic::postProcess() {
 if(domain->solInfo().statevectPodRom)
  stateImpl_->postProcess();
 if(domain->solInfo().residvectPodRom)
  resImpl_->postProcess();
 if(domain->solInfo().jacobvectPodRom)
  jacImpl_->postProcess();
}

void
SnapshotNonLinDynamic::saveMidTime(double t) {
 if(domain->solInfo().statevectPodRom)
        stateImpl_->lastMidTimeIs(t);
}

void
SnapshotNonLinDynamic::saveDelta(double dt) {
 if(domain->solInfo().statevectPodRom)
        stateImpl_->lastDeltaIs(dt);
}

void
SnapshotNonLinDynamic::saveStateSnapshot(const GeomState &state) {
 if(domain->solInfo().statevectPodRom)
        stateImpl_->stateSnapshotAdd(state);
}

void
SnapshotNonLinDynamic::handleResidualSnapshot(const Vector &snap) {
 if(domain->solInfo().residvectPodRom)
	resImpl_->handleResidualSnapshot(snap);
}

int
SnapshotNonLinDynamic::checkConvergence(int iteration, double normRes, Vector &residual, Vector &dv, double time) {
 if(domain->solInfo().jacobvectPodRom)
    jacImpl_->handleJacobianSnapshot();

  // Forward to hidden base class function
  return NonLinDynamic::checkConvergence(iteration, normRes, residual, dv, time);
}

} /* end namespace Rom */
