#ifndef ROM_VECBASISOPS_H
#define ROM_VECBASISOPS_H

#include "VecBasis.h"
#include "CholeskyUtils.h"

#include <Math.d/Vector.h>
#include <Math.d/SparseMatrix.h>
#include <Math.d/FullSquareMatrix.h>

#include <algorithm>
#include <memory>
#include <cassert>


namespace Rom {

// Returns the matrix-basis product
// (basis and result must refer to different objects)
template <typename Scalar>
const GenVecBasis<Scalar, GenVector> &
mult(const GenSparseMatrix<Scalar> &matrix, const GenVecBasis<Scalar, GenVector> &basis, GenVecBasis<Scalar, GenVector> &result) {
  assert(&basis != &result);

  typedef GenVecBasis<Scalar, GenVector> BasisType;
  typedef typename BasisType::iterator VecIt;

  result.dimensionIs(basis.vectorCount(), basis.vectorInfo());
  
  for (VecIt it = const_cast<BasisType &>(basis).begin(),
             it_end = const_cast<BasisType &>(basis).end(),
             jt = result.begin();
       it != it_end;
       ++it, ++jt) {
    const_cast<GenSparseMatrix<Scalar> &>(matrix).mult(*it, *jt);
  }

  return result;
}


// Returns the renormalized basis Phi with respect to the metric M (assumed symmetric positive semidefinite)
// (Phi, M) -> Phi * R^{-T} where (Phi^T * M * Phi) = R * R^T
// Distributed consistency: 1) requirements: Phi fully consistent, M does not assemble 2) guarantees: result fully consistent
template <typename Scalar>
const GenVecBasis<Scalar, GenVector> &
renormalized_basis(const GenSparseMatrix<Scalar> &metric, const GenVecBasis<Scalar, GenVector> &basis, GenVecBasis<Scalar, GenVector> &result) {
  // result <- M * Phi
  mult(metric, basis, result);

  // Build the lower-triangular part of the normal matrix
  // normalMatrix <- lower( (M * Phi)^T * Phi )
  const int vecCount = result.vectorCount();
  GenFullSquareMatrix<Scalar> normalMatrix(vecCount);

  for (int row = 0; row < vecCount; ++row) {
    const GenVector<Scalar> &dual = result[row];
    for (int col = 0; col <= row; ++col) {
      // sum-consistent * fully consistent
      normalMatrix[row][col] = dual * basis[col];
    }
  }

  cholesky_factor_lower(normalMatrix); // normalMatrix <- R where A = R * R^T
  inverse_triangular_lower(normalMatrix); // normalMatrix <- R^{-1}

  // result <- Phi * R^{-T}
  for (int row = 0; row < vecCount; ++row) {
    GenVector<Scalar> &target = result[row];
    int col = 0;
    target.linC(basis[col], normalMatrix[row][col]);
    for (col = 1; col <= row; ++col) {
      target.linAdd(normalMatrix[row][col], basis[col]);
    }
  }

  return result;
}

} // end namespace Rom

#endif /* ROM_VECBASISOPS_H */