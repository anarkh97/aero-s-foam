#include "MappedNodeDof6Buffer.h"

#include "RenumberingUtils.h"

#include <iterator>
#include <cstddef>

namespace Rom {

const double *
MappedNodeDof6Buffer::operator[](int iNode) const {
  std::map<int, int>::const_iterator it = underlyingNodeIndices_.find(iNode);
  return (it != underlyingNodeIndices_.end()) ? buffer_[it->second] : NULL;
}
  
void
MappedNodeDof6Buffer::initialize() {
  // Compute the inverse mapping (usual to underlying) 
  inverse_numbering(nodeIndices_.begin(), nodeIndices_.end(),
                    std::inserter(underlyingNodeIndices_, underlyingNodeIndices_.end()));
  
  // Resize internal buffer (underlying indexing)
  buffer_.sizeIs(nodeIndices_.size());
}

} /* end namespace Rom */
