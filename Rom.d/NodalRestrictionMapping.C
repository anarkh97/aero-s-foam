#include "NodalRestrictionMapping.h"

#include <Utils.d/dofset.h>

namespace Rom {

NodalRestrictionMapping::InfoType
NodalRestrictionMapping::extractOriginalInfo(const DofSetArray &dsa) {
  return const_cast<DofSetArray &>(dsa).size();
}

void
NodalRestrictionMapping::addSampleNode(int iNode, const DofSetArray &dsa) {
  static const int DOF_ID[] = { DofSet::Xdisp, DofSet::Ydisp, DofSet::Zdisp,
                                DofSet::Xrot,  DofSet::Yrot,  DofSet::Zrot  };
  
  for (const int *itDof = DOF_ID; itDof != DOF_ID + 6; ++itDof) {
    const int originLoc = const_cast<DofSetArray &>(dsa).locate(iNode, *itDof);
    if (originLoc >= 0) {
      originIndex_.push_back(originLoc);
    }
  }
}

std::ostream &
operator<<(std::ostream &out, const NodalRestrictionMapping &source) {
  out << "NodalRestrictionMapping: " << source.originInfo() << "->" << source.restrictedInfo() << " :";
  int index = 0;
  typedef std::vector<NodalRestrictionMapping::IndexType>::const_iterator idx_it;
  for (idx_it it = source.originIndex_.begin(); it != source.originIndex_.end(); ++it) {
    out << " " << *it << "->" << index++;
  }

  return out;
}

} /* end namespace Rom */
