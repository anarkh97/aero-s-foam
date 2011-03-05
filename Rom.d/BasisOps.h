#ifndef ROM_BASISOPS_H
#define ROM_BASISOPS_H

#include <cassert>

namespace Rom {

template <typename VecSetType, typename VecType>
inline
void
checkDimensions(const VecSetType &basis, const VecType &xFull, const VecType &xReduced) {
  assert(basis.size() == xFull.size());
  assert(basis.numVec() == xReduced.size());
}

template <typename VecSetType, typename VecType>
const VecType &
reduce(const VecSetType &basis, const VecType &xFull, VecType &xReduced) {
  checkDimensions(basis, xFull, xReduced);

  for (int i = 0; i < basis.numVec(); ++i) {
    xReduced[i] = basis[i] * xFull;
  }

  return xReduced;
}

template <typename VecSetType, typename VecType>
const VecType &
expand(const VecSetType &basis, const VecType &xReduced, VecType &xFull) {
  checkDimensions(basis, xFull, xReduced);

  xFull.zero();
  for (int i = 0; i < basis.numVec(); ++i) {
    xFull.linAdd(xReduced[i], basis[i]);
  }

  return xFull; 
}

} /* end namespace Rom */

#endif /* ROM_BASISOPS_H */
