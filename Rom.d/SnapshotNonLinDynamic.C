#include "SnapshotNonLinDynamic.h"

#include "SvdOrthogonalization.h"
#include "BasisOutputFile.h"

#include <Driver.d/Domain.h>

#include <Utils.d/dofset.h>
#include <Corotational.d/utilities.h>

#include <deque>

#include <cstddef>
#include <algorithm>

// Dummy class holding the implementation of SnapshotNonLinDynamic
struct SnapshotNonLinDynamicDetail : private SnapshotNonLinDynamic {
  class RawImpl : public Impl {
  public:
    virtual void stateSnapshotAdd(const GeomState &);
    virtual void residualSnapshotAdd(const GenVector<double> &);
    virtual void postProcess();

    int nodeCount() const { return domain_->numGlobalNodes(); }
    int vectorSize() const { return domain_->numUncon(); }

    virtual ~RawImpl();

    explicit RawImpl(Domain *);

  protected:
    void fillSnapBuffer(const double *origin);
  
    typedef double (*NodeBuffer)[6];
    const NodeBuffer snapBuffer() const { return snapBuffer_; }

  private:
    Domain * domain_;
    
    NodeBuffer snapBuffer_;
    int (*dofLocation_)[6];

  protected:
    BasisOutputFile stateSnapFile_;
    BasisOutputFile residualSnapFile_;
  };

  class SvdImpl : public RawImpl {
  public:
    virtual void stateSnapshotAdd(const GeomState &);
    virtual void residualSnapshotAdd(const GenVector<double> &);
    virtual void postProcess();

    explicit SvdImpl(Domain *);

    virtual ~SvdImpl();

  private:
    void orthoAndSave(const std::deque<Vector> &, BasisOutputFile &);

    std::deque<Vector> stateSnapshot_;
    std::deque<Vector> residualSnapshot_;

    const GeomState * refGeomState_;
    Vector increment_;

    SvdOrthogonalization svdSolver_;
  };

private:
  // Dummy constructor to avoid compilation failures
  SnapshotNonLinDynamicDetail(Domain *d) :
    SnapshotNonLinDynamic(d)
  {}
};

SnapshotNonLinDynamicDetail::RawImpl::RawImpl(Domain * domain) :
  domain_(domain),
  snapBuffer_(new double[nodeCount()][6]),
  dofLocation_(new int[nodeCount()][6]),
  stateSnapFile_("StateSnap", nodeCount()),      // TODO name
  residualSnapFile_("ResidualSnap", nodeCount()) // TODO name
{
  // Store location of each degree of freedom
  DofSetArray &cdsa = *(domain_->getCDSA());
  for (int iNode = 0; iNode < nodeCount(); ++iNode) {
    dofLocation_[iNode][0] = cdsa.locate(iNode, DofSet::Xdisp);
    dofLocation_[iNode][1] = cdsa.locate(iNode, DofSet::Ydisp);
    dofLocation_[iNode][2] = cdsa.locate(iNode, DofSet::Zdisp);
    dofLocation_[iNode][3] = cdsa.locate(iNode, DofSet::Xrot);
    dofLocation_[iNode][4] = cdsa.locate(iNode, DofSet::Yrot);
    dofLocation_[iNode][5] = cdsa.locate(iNode, DofSet::Zrot);
  }
}

SnapshotNonLinDynamicDetail::RawImpl::~RawImpl() {
  delete[] dofLocation_;
  delete[] snapBuffer_;
}

void
SnapshotNonLinDynamicDetail::RawImpl::postProcess() {
  stateSnapFile_.updateStateCountStatus();
  residualSnapFile_.updateStateCountStatus();
}

void
SnapshotNonLinDynamicDetail::RawImpl::stateSnapshotAdd(const GeomState &snap) {
  const CoordSet &refCoords = domain_->getNodes();

  for (int iNode = 0; iNode < nodeCount(); ++iNode) {
    // Translational dofs
    snapBuffer_[iNode][0] = snap[iNode].x - refCoords[iNode]->x;
    snapBuffer_[iNode][1] = snap[iNode].y - refCoords[iNode]->y;
    snapBuffer_[iNode][2] = snap[iNode].z - refCoords[iNode]->z;
    
    // Rotational dofs
    double rot[3];
    mat_to_vec(const_cast<double (*)[3]>(snap[iNode].R), rot);
    std::copy(rot, rot + 3, snapBuffer_[iNode] + 3);
  }
 
  stateSnapFile_.stateAdd(snapBuffer_);
}

void
SnapshotNonLinDynamicDetail::RawImpl::residualSnapshotAdd(const Vector &snap) {
  fillSnapBuffer(snap.data());
  residualSnapFile_.stateAdd(snapBuffer_);
}

void
SnapshotNonLinDynamicDetail::RawImpl::fillSnapBuffer(const double *snap) {
  for (int iNode = 0; iNode < nodeCount(); ++iNode) {
    for (int iDof = 0; iDof < 6; ++iDof) {
      const int loc = dofLocation_[iNode][iDof];
      snapBuffer_[iNode][iDof] = (loc >= 0) ? snap[loc] : 0.0;
    }
  }
}

SnapshotNonLinDynamicDetail::SvdImpl::SvdImpl(Domain * domain) :
  RawImpl(domain),
  refGeomState_(new GeomState(*domain->getDSA(), *domain->getCDSA(), domain->getNodes())),
  increment_(vectorSize())
{}

SnapshotNonLinDynamicDetail::SvdImpl::~SvdImpl() {
  delete refGeomState_;
}

void
SnapshotNonLinDynamicDetail::SvdImpl::stateSnapshotAdd(const GeomState &snap) {
  const_cast<GeomState &>(snap).diff(*refGeomState_, increment_);
  stateSnapshot_.push_back(increment_);
}

void
SnapshotNonLinDynamicDetail::SvdImpl::residualSnapshotAdd(const GenVector<double> &snap) {
  residualSnapshot_.push_back(snap);
}

void
SnapshotNonLinDynamicDetail::SvdImpl::postProcess() {
  orthoAndSave(stateSnapshot_, this->stateSnapFile_);
  orthoAndSave(residualSnapshot_, this->residualSnapFile_);

  RawImpl::postProcess();
}

void
SnapshotNonLinDynamicDetail::SvdImpl::orthoAndSave(const std::deque<Vector> & snapshots, BasisOutputFile & out) {
  svdSolver_.matrixSizeIs(vectorSize(), snapshots.size());
  
  int col = 0;
  const std::deque<Vector>::const_iterator itEnd = snapshots.end();
  for (std::deque<Vector>::const_iterator it = snapshots.begin(); it != itEnd; ++it) {
    const Vector &v = *it;
    std::copy(v.data(), v.data() + v.size(), svdSolver_.matrixCol(col));
    ++col;
  }

  svdSolver_.solve();
  
  const int stateCount = svdSolver_.singularValueCount();
  for (int iState = 0; iState < stateCount; ++iState) {
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

