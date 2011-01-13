#include "SnapshotNonLinDynamic.h"

#include "SvdOrthogonalization.h"
#include "BasisOutputFile.h"

#include <Driver.d/Domain.h>

#include <Utils.d/dofset.h>
#include <Corotational.d/utilities.h>

#include <list>

#include <cstddef>
#include <algorithm>

namespace {

class RawImpl : public SnapshotNonLinDynamic::Impl {
public:
  virtual void stateSnapshotAdd(const GeomState &);
  virtual void residualSnapshotAdd(const GenVector<double> &);
  virtual void postProcess();

  virtual ~RawImpl();

  explicit RawImpl(Domain *);

protected:
  void fillSnapBuffer(const double *origin);

  Domain * domain_;

  double (*snapBuffer_)[6];
  int (*dofLocation_)[6];

  BasisOutputFile stateSnapFile_;
  BasisOutputFile residualSnapFile_;
};

RawImpl::RawImpl(Domain * domain) :
  domain_(domain),
  snapBuffer_(new double[domain->numGlobalNodes()][6]),
  dofLocation_(new int[domain->numGlobalNodes()][6]),
  stateSnapFile_("StateSnap", domain->numGlobalNodes()),      // TODO name
  residualSnapFile_("ResidualSnap", domain->numGlobalNodes()) // TODO name
{
  // Store location of each degree of freedom
  const int nodeCount = domain_->numGlobalNodes();
  DofSetArray &cdsa = *(domain_->getCDSA());
  for (int iNode = 0; iNode < nodeCount; ++iNode) {
    dofLocation_[iNode][0] = cdsa.locate(iNode, DofSet::Xdisp);
    dofLocation_[iNode][1] = cdsa.locate(iNode, DofSet::Ydisp);
    dofLocation_[iNode][2] = cdsa.locate(iNode, DofSet::Zdisp);
    dofLocation_[iNode][3] = cdsa.locate(iNode, DofSet::Xrot);
    dofLocation_[iNode][4] = cdsa.locate(iNode, DofSet::Yrot);
    dofLocation_[iNode][5] = cdsa.locate(iNode, DofSet::Zrot);
  }
}

RawImpl::~RawImpl() {
  delete[] dofLocation_;
  delete[] snapBuffer_;
}

void
RawImpl::postProcess() {
  stateSnapFile_.updateStateCountStatus();
  residualSnapFile_.updateStateCountStatus();
}

void
RawImpl::stateSnapshotAdd(const GeomState &snap) {
  const CoordSet &refCoords = domain_->getNodes();
  const int nodeCount = domain_->numGlobalNodes();

  for (int iNode = 0; iNode < nodeCount; ++iNode) {
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
RawImpl::residualSnapshotAdd(const Vector &snap) {
  fillSnapBuffer(snap.data());
  residualSnapFile_.stateAdd(snapBuffer_);
}

void
RawImpl::fillSnapBuffer(const double *snap) {
  const int nodeCount = domain_->numGlobalNodes();
  
  for (int iNode = 0; iNode < nodeCount; ++iNode) {
    for (int iDof = 0; iDof < 6; ++iDof) {
      const int loc = dofLocation_[iNode][iDof];
      snapBuffer_[iNode][iDof] = (loc >= 0) ? snap[loc] : 0.0;
    }
  }
}

class SvdImpl : public RawImpl {
public:
  virtual void stateSnapshotAdd(const GeomState &);
  virtual void residualSnapshotAdd(const GenVector<double> &);
  virtual void postProcess();

  explicit SvdImpl(Domain *);

  virtual ~SvdImpl();

private:
  void orthoAndSave(const std::list<Vector> &, BasisOutputFile &);

  std::list<Vector> stateSnapshot_;
  std::list<Vector> residualSnapshot_;

  const GeomState * refGeomState_;
  Vector increment_;

  SvdOrthogonalization svdSolver_;
};

SvdImpl::SvdImpl(Domain * domain) :
  RawImpl(domain),
  refGeomState_(new GeomState(*domain->getDSA(), *domain->getCDSA(), domain->getNodes())),
  increment_(domain->numUncon())
{}

SvdImpl::~SvdImpl() {
  delete refGeomState_;
}

void
SvdImpl::stateSnapshotAdd(const GeomState &snap) {
  const_cast<GeomState &>(snap).diff(*refGeomState_, increment_);
  stateSnapshot_.push_back(increment_);
}

void
SvdImpl::residualSnapshotAdd(const GenVector<double> &snap) {
  residualSnapshot_.push_back(snap);
}

void
SvdImpl::postProcess() {
  orthoAndSave(stateSnapshot_, stateSnapFile_);
  orthoAndSave(residualSnapshot_, residualSnapFile_);

  RawImpl::postProcess();
}

void
SvdImpl::orthoAndSave(const std::list<Vector> & snapshots, BasisOutputFile & out) {
  svdSolver_.matrixSizeIs(domain->numUncon(), snapshots.size());
  
  int col = 0;
  for (std::list<Vector>::const_iterator it = snapshots.begin();
       it != snapshots.end(); ++it) {
    const Vector &v = *it;
    std::copy(v.data(), v.data() + v.size(), svdSolver_.matrixCol(col));
    ++col;
  }

  svdSolver_.solve();
  
  const int stateCount = svdSolver_.singularValueCount();
  const int nodeCount = domain_->numGlobalNodes();
 
  for (int iState = 0; iState < stateCount; ++iState) {
    fillSnapBuffer(svdSolver_.matrixCol(iState));
    out.stateAdd(snapBuffer_, svdSolver_.singularValue(iState));
  }
}

} // end anonymous namespace

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
    impl_.reset(new RawImpl(this->domain));
    break;
  case ORTHOGONAL:
    impl_.reset(new SvdImpl(this->domain));
    break;
  }
}

