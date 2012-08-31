#ifndef ROM_SPARSENONNEGATIVELEASTSQUARESSOLVER_H
#define ROM_SPARSENONNEGATIVELEASTSQUARESSOLVER_H

#include "SimpleBuffer.h"
#ifdef USE_STXXL
#include <stxxl/vector>
#else
#include <vector>
#endif

namespace Rom {

class SparseNonNegativeLeastSquaresSolver {
public:
  typedef double Scalar;
#ifdef USE_STXXL
  typedef stxxl::vector<Scalar> MatrixBufferType;
#else
  typedef std::vector<Scalar> MatrixBufferType;
#endif

  // Problem size
  int equationCount() const { return equationCount_; }
  int unknownCount() const { return unknownCount_; }
  
  void problemSizeIs(int eqnCount, int unkCount);

  double relativeTolerance() const { return relativeTolerance_; }
  void relativeToleranceIs(double relTol) { relativeTolerance_ = relTol; }

  // Buffers: Internal column-major ordering, zero-based indexing
  // Matrix buffer: [equationCount by unknownCount]
  Scalar matrixEntry(int row, int col) const;
  MatrixBufferType::const_iterator matrixColBuffer(int col) const;
  MatrixBufferType::const_iterator matrixBuffer() const;

  Scalar & matrixEntry(int row, int col);
  MatrixBufferType::iterator matrixColBuffer(int col);
  MatrixBufferType::iterator matrixBuffer();

  // Rhs buffer: [equationCount]
  Scalar rhsEntry(int row) const;
  const Scalar * rhsBuffer() const;
  
  Scalar & rhsEntry(int row);
  Scalar * rhsBuffer();

  // Primal solution buffer: [unknownCount]
  Scalar solutionEntry(int row) const;
  const Scalar * solutionBuffer() const;

  // Dual solution buffer: [unknownCount]
  Scalar dualSolutionEntry(int row) const;
  const Scalar * dualSolutionBuffer() const;

  // Error magnitude
  double errorMagnitude() const { return errorMagnitude_; }

  // Perform solve
  void solve();

  // Constructor
  SparseNonNegativeLeastSquaresSolver();

private:
  int equationCount_;
  int unknownCount_;
  int matrixLeadDim_;

  double relativeTolerance_;
  MatrixBufferType matrixBuffer_;

  SimpleBuffer<Scalar> rhsBuffer_;
  SimpleBuffer<Scalar> solutionBuffer_;
  SimpleBuffer<Scalar> dualSolutionBuffer_;

  double errorMagnitude_;

  // Disallow copy & assignment
  SparseNonNegativeLeastSquaresSolver(const SparseNonNegativeLeastSquaresSolver &);
  SparseNonNegativeLeastSquaresSolver &operator=(const SparseNonNegativeLeastSquaresSolver &);
};

// Inline member functions
inline
SparseNonNegativeLeastSquaresSolver::MatrixBufferType::const_iterator
SparseNonNegativeLeastSquaresSolver::matrixBuffer() const {
  return matrixBuffer_.begin();
} 

inline
SparseNonNegativeLeastSquaresSolver::MatrixBufferType::const_iterator
SparseNonNegativeLeastSquaresSolver::matrixColBuffer(int col) const {
  return matrixBuffer() + (col * matrixLeadDim_);
}

inline
SparseNonNegativeLeastSquaresSolver::Scalar
SparseNonNegativeLeastSquaresSolver::matrixEntry(int row, int col) const {
  return matrixColBuffer(col)[row];
}

inline
SparseNonNegativeLeastSquaresSolver::MatrixBufferType::iterator
SparseNonNegativeLeastSquaresSolver::matrixBuffer() {
  return matrixBuffer_.begin();
} 

inline
SparseNonNegativeLeastSquaresSolver::MatrixBufferType::iterator
SparseNonNegativeLeastSquaresSolver::matrixColBuffer(int col) {
  return matrixBuffer() + (col * matrixLeadDim_);
}

inline
SparseNonNegativeLeastSquaresSolver::Scalar &
SparseNonNegativeLeastSquaresSolver::matrixEntry(int row, int col) {
  return matrixColBuffer(col)[row];
}

inline
const SparseNonNegativeLeastSquaresSolver::Scalar *
SparseNonNegativeLeastSquaresSolver::rhsBuffer() const {
  return rhsBuffer_.array();
} 

inline
SparseNonNegativeLeastSquaresSolver::Scalar
SparseNonNegativeLeastSquaresSolver::rhsEntry(int row) const {
  return rhsBuffer()[row];
}

inline
SparseNonNegativeLeastSquaresSolver::Scalar *
SparseNonNegativeLeastSquaresSolver::rhsBuffer() {
  return rhsBuffer_.array();
} 

inline
SparseNonNegativeLeastSquaresSolver::Scalar &
SparseNonNegativeLeastSquaresSolver::rhsEntry(int row) {
  return rhsBuffer()[row];
}

inline
const SparseNonNegativeLeastSquaresSolver::Scalar *
SparseNonNegativeLeastSquaresSolver::solutionBuffer() const {
  return solutionBuffer_.array();
} 

inline
SparseNonNegativeLeastSquaresSolver::Scalar
SparseNonNegativeLeastSquaresSolver::solutionEntry(int row) const {
  return solutionBuffer()[row];
}

inline
const SparseNonNegativeLeastSquaresSolver::Scalar *
SparseNonNegativeLeastSquaresSolver::dualSolutionBuffer() const {
  return dualSolutionBuffer_.array();
} 

inline
SparseNonNegativeLeastSquaresSolver::Scalar
SparseNonNegativeLeastSquaresSolver::dualSolutionEntry(int row) const {
  return dualSolutionBuffer()[row];
}

} // end namespace Rom

#endif /* ROM_SPARSENONNEGATIVELEASTSQUARESSOLVER_H */
