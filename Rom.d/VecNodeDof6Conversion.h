#ifndef ROM_VECNODEDOF6CONVERSION_H
#define ROM_VECNODEDOF6CONVERSION_H

#include "SimpleBuffer.h"

#include <cassert>

class DofSetArray;

namespace Rom {

class VecNodeDof6Conversion {
public:
  int nodeCount() const { return nodeCount_; }
  int vectorSize() const { return vectorSize_; }

  template <typename NodeDofs6Type, typename VecType>
  const NodeDofs6Type &nodeDof6(const VecType &origin, NodeDofs6Type &target) const;

  template <typename NodeDofs6Type, typename VecType>
  const VecType &vector(const NodeDofs6Type &origin, VecType &target) const;

  explicit VecNodeDof6Conversion(const DofSetArray &);

private:
  int nodeCount_;
  int vectorSize_;

  typedef SimpleBuffer<int[6]> DofLocation;
  DofLocation dofLocation_;

  // Disallow copy and assignment
  VecNodeDof6Conversion(const VecNodeDof6Conversion &);
  VecNodeDof6Conversion &operator=(const VecNodeDof6Conversion &);
};

template <typename NodeDofs6Type, typename VecType>
const NodeDofs6Type &
VecNodeDof6Conversion::nodeDof6(const VecType &origin, NodeDofs6Type &target) const {
  for (int iNode = 0; iNode < nodeCount(); ++iNode) {
    for (int iDof = 0; iDof < 6; ++iDof) {
      const int loc = dofLocation_[iNode][iDof];
      target[iNode][iDof] = (loc >= 0) ? origin[loc] : 0.0;
    }
  }

  return target;
}

template <typename NodeDofs6Type, typename VecType>
const VecType &
VecNodeDof6Conversion::vector(const NodeDofs6Type &origin, VecType &target) const {
  for (int iNode = 0; iNode < nodeCount(); ++iNode) {
    for (int iDof = 0; iDof < 6; ++iDof) {
      const int loc = dofLocation_[iNode][iDof];
      if (loc >= 0) {
        target[loc] = origin[iNode][iDof];
      }
    }
  }

  return target;
}

} /* end namespace Rom */

#endif /* ROM_VECNODEDOF6CONVERSION_H */
