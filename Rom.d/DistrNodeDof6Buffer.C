#include "DistrNodeDof6Buffer.h"

#include "RenumberingUtils.h"

#include <algorithm>
#include <iterator>
#include <cstddef>

namespace Rom {

const double *
DistrNodeDof6Buffer::operator[](int globalNodeIdx) const {
  std::map<int, int>::const_iterator it = localNodeIndices_.find(globalNodeIdx);
  return (it != localNodeIndices_.end()) ? buffer_[it->second] : NULL;
}
  
void
DistrNodeDof6Buffer::initialize() {
  // Sort and make sure every index is unique (local to global)
  std::sort(globalNodeIndices_.begin(), globalNodeIndices_.end());
  globalNodeIndices_.erase(std::unique(globalNodeIndices_.begin(), globalNodeIndices_.end()), globalNodeIndices_.end());

  // Compute the inverse mapping (global to local) 
  inverse_numbering(globalNodeIndices_.begin(), globalNodeIndices_.end(),
                    std::inserter(localNodeIndices_, localNodeIndices_.end()));
  
  // Resize internal buffer (local indexing)
  buffer_.sizeIs(localNodeCount());
}

} /* end namespace Rom */
