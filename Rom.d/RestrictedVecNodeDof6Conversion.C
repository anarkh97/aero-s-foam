#include "RestrictedVecNodeDof6Conversion.h"

namespace Rom {

void
RestrictedVecNodeDof6Conversion::initialize(const DofSetArray &dsa) {
  for (int iNode = 0; iNode < nodeCount(); ++iNode) {
    for (int iDof = 0; iDof < DOF_ID_COUNT; ++iDof) {
      const NodeDof::DofType dofId = DOF_ID[iDof];
      const int loc = const_cast<DofSetArray &>(dsa).locate(iNode, dofId);

      dofLocation_[iNode][iDof] = loc;
      if (loc >= 0) {
        locationId_[loc] = NodeDof(iNode, dofId);
      }
    }
  }
}

} /* end namespace Rom */
