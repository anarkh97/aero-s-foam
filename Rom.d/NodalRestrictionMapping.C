#include "NodalRestrictionMapping.h"

#include <Utils.d/dofset.h>

const int
NodalRestrictionMapping::DOF_ID[6] = { DofSet::Xdisp, DofSet::Ydisp, DofSet::Zdisp,
                                       DofSet::Xrot,  DofSet::Yrot,  DofSet::Zrot  };

NodalRestrictionMapping::InfoType
NodalRestrictionMapping::extractOriginalInfo(const DofSetArray &dsa) {
  return const_cast<DofSetArray &>(dsa).size();
}

void
NodalRestrictionMapping::addSampleNode(int iNode, const DofSetArray &dsa) {
  for (const int *itDof = DOF_ID; itDof != DOF_ID + 6; ++itDof) {
    const int originLoc = const_cast<DofSetArray &>(dsa).locate(iNode, *itDof);
    if (originLoc >= 0) {
      originIndex_.push_back(originLoc);
    }
  }
}

