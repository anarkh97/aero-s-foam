#include "NodalRestrictionMapping.h"
#include "DofSetUtils.h"

#include <Utils.d/dofset.h>

namespace Rom {

NodalRestrictionMapping::InfoType
NodalRestrictionMapping::extractOriginalInfo(const DofSetArray &dsa) {
  return const_cast<DofSetArray &>(dsa).size();
}

void
NodalRestrictionMapping::addSampleNode(int iNode, const DofSetArray &dsa) {
  for (const NodeDof::DofType *itDof = DOF_ID, *itDofEnd = DOF_ID + DOF_ID_COUNT; itDof != itDofEnd; ++itDof) {
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
