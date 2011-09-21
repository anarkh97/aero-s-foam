#include "VecNodeDof6Conversion.h"
#include "DofSetUtils.h"

#include <Utils.d/dofset.h>

namespace Rom {

VecNodeDof6Conversion::VecNodeDof6Conversion(const DofSetArray &dsa) :
  nodeCount_(const_cast<DofSetArray &>(dsa).numNodes()),
  vectorSize_(const_cast<DofSetArray &>(dsa).size()),
  dofLocation_(nodeCount())
{
  for (int iNode = 0, iNodeEnd = nodeCount(); iNode < iNodeEnd; ++iNode) {
    for (int iDof = 0; iDof != DOF_ID_COUNT; ++iDof) {
      const int loc = const_cast<DofSetArray &>(dsa).locate(iNode, DOF_ID[iDof]);
      dofLocation_[iNode][iDof] = loc;
    }
  }
}

} /* end namespace Rom */
