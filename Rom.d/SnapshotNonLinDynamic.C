#include "SnapshotNonLinDynamic.h"

#include "SvdOrthogonalization.h"
#include "BasisOutputFile.h"

#include <Driver.d/Domain.h>

#include <Utils.d/dofset.h>
#include <Corotational.d/utilities.h>

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

private:
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
  const int nodeCount = domain_->numGlobalNodes();
  
  for (int iNode = 0; iNode < nodeCount; ++iNode) {
    for (int iDof = 0; iDof < 6; ++iDof) {
      const int loc = dofLocation_[iNode][iDof];
      snapBuffer_[iNode][iDof] = (loc >= 0) ? snap[loc] : 0.0;
    }
  }

  residualSnapFile_.stateAdd(snapBuffer_);
}

} // end anonymous namespace

SnapshotNonLinDynamic::SnapshotNonLinDynamic(Domain *domain, BasisType outputBasisType) :
  NonLinDynamic(domain),
  outputBasisType_(outputBasisType),
  impl_(NULL)
{}

void
SnapshotNonLinDynamic::preProcess() {
  NonLinDynamic::preProcess();

  // TODO: add case outputBasisType_ == ORTHOGONAL
  impl_.reset(new RawImpl(this->domain));
}

