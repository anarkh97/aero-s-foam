#include "QrPseudoInversion.h"

#include <Math.d/Vector.h>
#include <Utils.d/linkfc.h>

#include <algorithm>
#include <stdexcept>
#include <cassert>

extern "C" {
  // Perform the QR factorization
  void _FORTRAN(dgeqrf)(const int *m, const int *n, double *a, const int *lda,
                        double *tau, double *work, const int *lwork, int *info);

  // Solves trinagular system
  void _FORTRAN(dtrsm)(const char *side, const char *uplo, const char *transa, const char* diag,
                       const int *m, const int *n, const double *alpha, const double *a, const int *lda,
                       double *b, const int *ldb);
}

namespace Rom {

// Computes (Phi^{+})^T = (Phi * R^{-1}) * R^{-T}
const VecBasis &
QrPseudoInversion::operator()(VecBasis &basis) const {
  const int vectorSize = basis.vectorSize();
  const int vectorCount = basis.vectorCount();

  if (vectorSize < vectorCount) {
    throw std::domain_error("Matrix must have full column rank");
  }

  const int matrixBufferSize = vectorCount * vectorSize;
  double * basisBufferAddr = &basis[0][0];

  ScalarBuffer matrixBuffer(matrixBufferSize);
  ScalarBuffer tauBuffer(vectorCount);

  std::copy(basisBufferAddr, basisBufferAddr + matrixBufferSize, matrixBuffer.array());

  int info;

  // Computes Phi = Q * R
  const int lworkQuery = -1;
  double lworkAns;
  _FORTRAN(dgeqrf)(&vectorSize, &vectorCount, matrixBuffer.array(), &vectorSize,
                   tauBuffer.array(), &lworkAns, &lworkQuery, &info);
  assert(info == 0);

  const int lwork = static_cast<int>(lworkAns);
  ScalarBuffer workBuffer(lwork);
  _FORTRAN(dgeqrf)(&vectorSize, &vectorCount, matrixBuffer.array(), &vectorSize,
                   tauBuffer.array(), workBuffer.array(), &lwork, &info);
  assert(info == 0);

  const char right = 'R';
  const char upper = 'U';
  const char lower = 'L';
  const char trans = 'T';
  const char no    = 'N'; // Not transposed, non-unit
  const double one = 1.0;

  // Computes Phi * R^{-1} 
  _FORTRAN(dtrsm)(&right, &upper, &no, &no,
                  &vectorSize, &vectorCount, &one, matrixBuffer.array(), &vectorSize,
                  basisBufferAddr, &vectorSize);

  // Computes (Phi * R^{-1}) * R^{-T} 
  _FORTRAN(dtrsm)(&right, &upper, &trans, &no,
                  &vectorSize, &vectorCount, &one, matrixBuffer.array(), &vectorSize,
                  basisBufferAddr, &vectorSize);
  
  return basis;
}

} /* end namespace Rom */
