#ifndef ROM_BASISOPS_H
#define ROM_BASISOPS_H

#include <Math.d/Vector.h>
#include <Math.d/VectorSet.h>

#include <cassert>

template <typename Scalar>
inline
void
checkDimensions(const GenVectorSet<Scalar> &basis, const GenVector<Scalar> &xFull, const GenVector<Scalar> &xReduced) {
  assert(basis.size() == xFull.size());
  assert(basis.numVec() == xReduced.size());
}

template <typename Scalar>
const GenVector<Scalar> &
reduce(const GenVectorSet<Scalar> &basis, const GenVector<Scalar> &xFull, GenVector<Scalar> &xReduced) {
  checkDimensions(basis, xFull, xReduced);

  for (int i = 0; i < basis.size(); ++i) {
    xReduced[i] = basis[i] * xFull;
  }

  return xReduced;
}

template <typename Scalar>
const GenVector<Scalar> &
expand(const GenVectorSet<Scalar> &basis, const GenVector<Scalar> &xReduced, GenVector<Scalar> &xFull) {
  checkDimensions(basis, xFull, xReduced);

  xFull.zero();
  for (int i = 0; i < basis.size(); ++i) {
    xFull.linAdd(xReduced[i], basis[i]);
  }

  return xFull; 
}

#endif /* ROM_BASISOPS_H */
