#include "SnapshotNonLinDynamic.h"

#include "SvdOrthogonalization.h"
#include "VecNodeDof6Conversion.h"
#include "BasisBinaryFile.h"
#include "NodeDof6Buffer.h"
#include "FileNameInfo.h"

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
  class RawImpl : public Impl {
  public:
    virtual void stateSnapshotAdd(const GeomState &);
    virtual void postProcess();

    int dofSetNodeCount() const { return converter_.dofSetNodeCount(); }
    int vectorSize() const { return converter_.vectorSize(); }

    explicit RawImpl(Domain *, BasisId::Level level = BasisId::SNAPSHOTS);

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

  protected:
    FileNameInfo fileInfo_;
    BasisBinaryOutputFile stateSnapFile_;
  };

  class SvdImpl : public RawImpl {
  public:
    virtual void stateSnapshotAdd(const GeomState &);
    virtual void postProcess();

    explicit SvdImpl(Domain *);

  private:
    void orthoAndSave(const std::deque<Vector> &, BasisBinaryOutputFile &);

    std::deque<Vector> stateSnapshot_;

    std::auto_ptr<const GeomState> refGeomState_;
    Vector increment_;

    SvdOrthogonalization svdSolver_;
  };

private:
  // Dummy constructor to avoid compilation failures
  SnapshotNonLinDynamicDetail(Domain *d) :
    SnapshotNonLinDynamic(d)
  {}
};

SnapshotNonLinDynamicDetail::RawImpl::RawImpl(Domain * domain, BasisId::Level level) :
  domain_(domain),
  converter_(*domain->getCDSA()),
  snapBuffer_(dofSetNodeCount()),
  fileInfo_(),
  stateSnapFile_(BasisFileId(fileInfo_, BasisId::STATE, level), dofSetNodeCount())
{}

void
SnapshotNonLinDynamicDetail::RawImpl::postProcess() {
  // Nothing to do
}

template <typename VecType>
inline
void
SnapshotNonLinDynamicDetail::RawImpl::fillSnapBuffer(const VecType &snap) {
  converter_.paddedNodeDof6(snap, snapBuffer_);
}

void
SnapshotNonLinDynamicDetail::RawImpl::stateSnapshotAdd(const GeomState &snap) {
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
 
  stateSnapFile_.stateAdd(snapBuffer_);
}

SnapshotNonLinDynamicDetail::SvdImpl::SvdImpl(Domain * domain) :
  RawImpl(domain, BasisId::POD),
  refGeomState_(new GeomState(*domain->getDSA(), *domain->getCDSA(), domain->getNodes())),
  increment_(vectorSize())
{}

void
SnapshotNonLinDynamicDetail::SvdImpl::stateSnapshotAdd(const GeomState &snap) {
  const_cast<GeomState &>(snap).diff(*refGeomState_, increment_);
  stateSnapshot_.push_back(increment_);
}

void
SnapshotNonLinDynamicDetail::SvdImpl::postProcess() {
  orthoAndSave(stateSnapshot_, this->stateSnapFile_);

  RawImpl::postProcess();
}

void
SnapshotNonLinDynamicDetail::SvdImpl::orthoAndSave(const std::deque<Vector> & snapshots, BasisBinaryOutputFile & out) {
  svdSolver_.matrixSizeIs(vectorSize(), snapshots.size());
  
  int col = 0;
  const std::deque<Vector>::const_iterator itEnd = snapshots.end();
  for (std::deque<Vector>::const_iterator it = snapshots.begin(); it != itEnd; ++it) {
    const Vector &v = *it;
    std::copy(v.data(), v.data() + v.size(), svdSolver_.matrixCol(col));
    ++col;
  }

  svdSolver_.solve();
  
  const int orthoBasisDim = maxSizePodRom() ?
                            std::min(maxSizePodRom(), svdSolver_.singularValueCount()) :
                            svdSolver_.singularValueCount();
  
  for (int iState = 0; iState < orthoBasisDim; ++iState) {
    fillSnapBuffer(svdSolver_.matrixCol(iState));
    out.stateAdd(snapBuffer(), svdSolver_.singularValue(iState));
  }
}

SnapshotNonLinDynamic::SnapshotNonLinDynamic(Domain *domain) :
  NonLinDynamic(domain),
  outputBasisType_(domain->solInfo().svdPodRom ? ORTHOGONAL : RAW),
  impl_(NULL)
{}

void
SnapshotNonLinDynamic::preProcess() {
  NonLinDynamic::preProcess();

  switch (outputBasisType_) {
  case RAW:
    impl_.reset(new SnapshotNonLinDynamicDetail::RawImpl(this->domain));
    break;
  case ORTHOGONAL:
    impl_.reset(new SnapshotNonLinDynamicDetail::SvdImpl(this->domain));
    break;
  }
}

} /* end namespace Rom */
