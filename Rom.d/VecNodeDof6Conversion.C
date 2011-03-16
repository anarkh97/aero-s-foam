#include "VecNodeDof6Conversion.h"

#include <Utils.d/dofset.h>

namespace Rom {

VecNodeDof6Conversion::VecNodeDof6Conversion(const DofSetArray &dsa) :
  nodeCount_(const_cast<DofSetArray &>(dsa).numNodes()),
  vectorSize_(const_cast<DofSetArray &>(dsa).size()),
  locationId_(vectorSize()),
  dofLocation_(new int[nodeCount()][6])
{
  static const int DOF_ID[] = { DofSet::Xdisp, DofSet::Ydisp, DofSet::Zdisp,
                                DofSet::Xrot,  DofSet::Yrot,  DofSet::Zrot  };

  for (int iNode = 0; iNode < nodeCount(); ++iNode) {
    for (int iDof = 0; iDof < 6; ++iDof) {
      const NodeDof::DofType dofId = DOF_ID[iDof];
      const int loc = const_cast<DofSetArray &>(dsa).locate(iNode, dofId);

      dofLocation_[iNode][iDof] = loc;
      if (loc >= 0) {
        locationId_[loc] = NodeDof(iNode, dofId);
      }
    }
  }
}

NodeDof
VecNodeDof6Conversion::nodeDof(int vecLoc) const {
  assert(vecLoc >= 0 && vecLoc < vectorSize());
  return locationId_[vecLoc]; 
} 

VecNodeDof6Conversion::~VecNodeDof6Conversion() {
  delete[] dofLocation_;
}

} /* end namespace Rom */
