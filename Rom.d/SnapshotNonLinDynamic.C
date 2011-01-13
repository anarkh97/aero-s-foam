#include <Rom.d/SnapshotNonLinDynamic.h>

#include <Driver.d/Domain.h>

#include <Utils.d/dofset.h>
#include <Corotational.d/utilities.h>

#include <cstddef>
#include <algorithm>

SnapshotNonLinDynamic::SnapshotNonLinDynamic(Domain *domain) :
  NonLinDynamic(domain),
  snapBuffer_(NULL),
  dofLocation_(NULL)
{}

SnapshotNonLinDynamic::~SnapshotNonLinDynamic() {
  delete[] dofLocation_;
  delete[] snapBuffer_;
}

void
SnapshotNonLinDynamic::preProcess() {
  NonLinDynamic::preProcess();
 
  const int nodeCount = this->domain->numGlobalNodes();

  snapBuffer_  = new double[nodeCount][6];
  dofLocation_ = new int[nodeCount][6];

  // Store location of each degree of freedom
  DofSetArray &cdsa = *(this->domain->getCDSA());
  for (int iNode = 0; iNode < nodeCount; ++iNode) {
    dofLocation_[iNode][0] = cdsa.locate(iNode, DofSet::Xdisp);
    dofLocation_[iNode][1] = cdsa.locate(iNode, DofSet::Ydisp);
    dofLocation_[iNode][2] = cdsa.locate(iNode, DofSet::Zdisp);
    dofLocation_[iNode][3] = cdsa.locate(iNode, DofSet::Xrot);
    dofLocation_[iNode][4] = cdsa.locate(iNode, DofSet::Yrot);
    dofLocation_[iNode][5] = cdsa.locate(iNode, DofSet::Zrot);
  }

  snapFile_[STATE_SNAP].reset(new BasisOutputFile("StateSnap", nodeCount)); // TODO name
  snapFile_[RESIDUAL_SNAP].reset(new BasisOutputFile("ResidualSnap", nodeCount)); // TODO name
}

void
SnapshotNonLinDynamic::saveStateSnapshot(const GeomState &snap) {
  const CoordSet &refCoords = this->domain->getNodes();
  const int nodeCount = this->domain->numGlobalNodes();

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
 
  snapFile_[STATE_SNAP]->stateAdd(snapBuffer_);
}

void
SnapshotNonLinDynamic::saveResidualSnapshot(const Vector &snap) {
  const int nodeCount = this->domain->numGlobalNodes();
  
  for (int iNode = 0; iNode < nodeCount; ++iNode) {
    for (int iDof = 0; iDof < 6; ++iDof) {
      const int loc = dofLocation_[iNode][iDof];
      snapBuffer_[iNode][iDof] = (loc >= 0) ? snap[loc] : 0.0;
    }
  }

  snapFile_[RESIDUAL_SNAP]->stateAdd(snapBuffer_); 
}
