#ifndef ROM_RESTRICTEDVECNODEDOF6CONVERSION_H
#define ROM_RESTRICTEDVECNODEDOF6CONVERSION_H

#include "DofSetUtils.h"
#include "SimpleBuffer.h"

#include <Utils.d/dofset.h>

#include <vector>
#include <cassert>

namespace Rom {

class RestrictedVecNodeDof6Conversion {
public:
  int nodeCount() const { return nodeCount_; }
  int vectorSize() const { return vectorSize_; }

  template <typename NodeDofs6Type, typename VecType>
  const NodeDofs6Type &paddedNodeDof6(const VecType &origin, NodeDofs6Type &target) const;

  template <typename NodeDofs6Type, typename VecType>
  const VecType &paddedVector(const NodeDofs6Type &origin, VecType &target) const;

  template <typename BoolFwdIt>
  RestrictedVecNodeDof6Conversion(const DofSetArray &dsa, BoolFwdIt masterFlagBegin,
                                                          BoolFwdIt masterFlagEnd);

private:
  int nodeCount_;
  int vectorSize_;

  std::vector<NodeDof> locationId_;
  
  typedef SimpleBuffer<int[6]> DofLocation;
  DofLocation dofLocation_;

  std::vector<bool> masterFlag;

  void initialize(const DofSetArray &dsa);

  // Disallow copy and assignment
  RestrictedVecNodeDof6Conversion(const RestrictedVecNodeDof6Conversion &);
  RestrictedVecNodeDof6Conversion &operator=(const RestrictedVecNodeDof6Conversion &);
};

template <typename NodeDofs6Type, typename VecType>
const NodeDofs6Type &
RestrictedVecNodeDof6Conversion::paddedNodeDof6(const VecType &origin, NodeDofs6Type &target) const {
  for (int iNode = 0; iNode < nodeCount(); ++iNode) {
    if (masterFlag[iNode]) {
      for (int iDof = 0; iDof < 6; ++iDof) {
        const int loc = dofLocation_[iNode][iDof];
        target[iNode][iDof] = (loc >= 0) ? origin[loc] : 0.0;
      }
    }
  }

  return target;
}

template <typename NodeDofs6Type, typename VecType>
const VecType &
RestrictedVecNodeDof6Conversion::paddedVector(const NodeDofs6Type &origin, VecType &target) const {
  for (int iNode = 0; iNode < nodeCount(); ++iNode) {
    if (masterFlag[iNode]) {
      for (int iDof = 0; iDof < 6; ++iDof) {
        const int loc = dofLocation_[iNode][iDof];
        if (loc >= 0) {
          target[loc] = origin[iNode][iDof];
        }
      }
    } else {
      for (int iDof = 0; iDof < 6; ++iDof) {
        const int loc = dofLocation_[iNode][iDof];
        if (loc >= 0) {
          target[loc] = 0.0;
        }
      }
    }
  }
  return target;
}

template <typename BoolFwdIt>
RestrictedVecNodeDof6Conversion::RestrictedVecNodeDof6Conversion(const DofSetArray &dsa,
                                                                 BoolFwdIt masterFlagBegin,
                                                                 BoolFwdIt masterFlagEnd) :
  nodeCount_(const_cast<DofSetArray &>(dsa).numNodes()),
  vectorSize_(const_cast<DofSetArray &>(dsa).size()),
  locationId_(vectorSize()),
  dofLocation_(nodeCount()),
  masterFlag(masterFlagBegin, masterFlagEnd)
{
  initialize(dsa);
}

} /* end namespace Rom */

#endif /* ROM_RESTRICTEDVECNODEDOF6CONVERSION_H */
