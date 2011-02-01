#ifndef ROM_LEASTSQUARESSOLVER_H
#define ROM_LEASTSQUARESSOLVER_H

#include "SimpleBuffer.h"

#include <stdexcept> 

template <typename Scalar>
class LeastSquaresSolver {
public:
  // Problem size
  int equationCount()    const { return equationCount_;    }
  int unknownCount()     const { return unknownCount_;     }
  int largestDimension() const { return largestDimension_; } // = max(eqnCount, unknCount)
  int rhsCount()         const { return rhsCount_;         }

  void problemSizeIs(int eqnCount, int unknownCount, int rhsCount = 1);
  
  // Buffers: Internal column-major ordering, zero-based indexing
  // Matrix buffer: [equationCount by unknownCount]
  Scalar matrixEntry(int row, int col) const;
  const Scalar * matrixColBuffer(int col) const;
  const Scalar * matrixBuffer() const;
  
  Scalar & matrixEntry(int row, int col);
  Scalar * matrixColBuffer(int col);
  Scalar * matrixBuffer();
  
  // Local Rhs/Solution buffer: [max(equationCount, unknownCount) by rhsCount]
  Scalar rhsEntry(int row) const; // First right-hand side
  Scalar rhsEntry(int rank, int row) const;
  const Scalar * rhsBuffer(int rank) const;
  const Scalar * rhsBuffer() const;
  
  Scalar & rhsEntry(int row); // First right-hand side
  Scalar & rhsEntry(int rank, int row);
  Scalar * rhsBuffer(int rank);
  Scalar * rhsBuffer();

  // Operations
  void factor();
  void solve();

  // Ctor 
  LeastSquaresSolver();

private:
  static void assertOverdeterminedSystem(int eqnCount, int unknownCount);

  int equationCount_, unknownCount_;
  int largestDimension_, rhsCount_;
  int matrixLeadDim_, rhsLeadDim_;

  typedef SimpleBuffer<Scalar> ScalarBuffer;
  
  ScalarBuffer matrixBuffer_;
  ScalarBuffer rhsBuffer_;
  ScalarBuffer tauBuffer_;

  // Disallow copy and assignment
  LeastSquaresSolver(const LeastSquaresSolver &);
  LeastSquaresSolver &operator=(const LeastSquaresSolver &);
};

// Inline member functions

template <typename Scalar>
inline
const Scalar *
LeastSquaresSolver<Scalar>::matrixBuffer() const {
  return matrixBuffer_.array();
} 

template <typename Scalar>
inline
const Scalar *
LeastSquaresSolver<Scalar>::matrixColBuffer(int col) const {
  return matrixBuffer() + (col * matrixLeadDim_);
}

template <typename Scalar>
inline
Scalar
LeastSquaresSolver<Scalar>::matrixEntry(int row, int col) const {
  return matrixColBuffer(col)[row];
}

template <typename Scalar>
inline
Scalar *
LeastSquaresSolver<Scalar>::matrixBuffer() {
  return matrixBuffer_.array();
} 

template <typename Scalar>
inline
Scalar *
LeastSquaresSolver<Scalar>::matrixColBuffer(int col) {
  return matrixBuffer() + (col * matrixLeadDim_);
}

template <typename Scalar>
inline
Scalar &
LeastSquaresSolver<Scalar>::matrixEntry(int row, int col) {
  return matrixColBuffer(col)[row];
}

template <typename Scalar>
inline
const Scalar *
LeastSquaresSolver<Scalar>::rhsBuffer() const {
  return rhsBuffer_.array();
} 

template <typename Scalar>
inline
const Scalar *
LeastSquaresSolver<Scalar>::rhsBuffer(int rank) const {
  return rhsBuffer() + (rank * rhsLeadDim_);
}

template <typename Scalar>
inline
Scalar
LeastSquaresSolver<Scalar>::rhsEntry(int rank, int row) const {
  return rhsBuffer(rank)[row];
}

template <typename Scalar>
inline
Scalar
LeastSquaresSolver<Scalar>::rhsEntry(int row) const {
  return rhsEntry(0, row);
}

template <typename Scalar>
inline
Scalar *
LeastSquaresSolver<Scalar>::rhsBuffer() {
  return rhsBuffer_.array();
} 

template <typename Scalar>
inline
Scalar *
LeastSquaresSolver<Scalar>::rhsBuffer(int rank) {
  return rhsBuffer() + (rank * rhsLeadDim_);
}

template <typename Scalar>
inline
Scalar &
LeastSquaresSolver<Scalar>::rhsEntry(int rank, int row) {
  return rhsBuffer(rank)[row];
}

template <typename Scalar>
inline
Scalar &
LeastSquaresSolver<Scalar>::rhsEntry(int row) {
  return rhsEntry(0, row);
}

// Implementation

template <typename Scalar>
LeastSquaresSolver<Scalar>::LeastSquaresSolver() :
  equationCount_(0),
  unknownCount_(0),
  largestDimension_(0),
  rhsCount_(0),
  matrixLeadDim_(1),
  rhsLeadDim_(1)
{}

template <typename Scalar>
void
LeastSquaresSolver<Scalar>::problemSizeIs(int eqnCount, int unknownCount, int rhsCount) {
  assertOverdeterminedSystem(eqnCount, unknownCount);

  equationCount_    = eqnCount;
  unknownCount_     = unknownCount;
  largestDimension_ = std::max(equationCount_, unknownCount_);
  rhsCount_         = rhsCount;

  matrixLeadDim_ = std::max(equationCount_, 1);
  rhsLeadDim_    = std::max(largestDimension_, 1);

  matrixBuffer_.sizeIs(equationCount_ * unknownCount_);
  rhsBuffer_.sizeIs(largestDimension_ * rhsCount_);
  tauBuffer_.sizeIs(unknownCount_);
}

template <typename Scalar>
void
LeastSquaresSolver<Scalar>::assertOverdeterminedSystem(int eqnCount, int unknownCount) {
  if (eqnCount < unknownCount) {
    throw std::domain_error("Underdetermined system");
  }
}

#endif /* ROM_LEASTSQUARESSOLVER_H */
