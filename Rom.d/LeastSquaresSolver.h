#ifndef ROM_LEASTSQUARESSOLVER_H
#define ROM_LEASTSQUARESSOLVER_H

#include "SimpleBuffer.h"

#include <stdexcept> 

namespace Rom {

namespace LeastSquares {
  // Describes the computational status of a LeastSquaresSolver
  enum Status { READY, FACTORED, PROJECTED, SOLVED };

  // READY     => Matrix & rhs buffers available for input
  // FACTORED  => Matrix buffer used for factorization, rhs buffer available for input
  // PROJECTED => Matrix buffer used for factorization, rhs buffer holds the projection Q^T * rhs
  // SOLVED    => Matrix buffer used for factorization, rhs buffer holds the solution R^{-1} * Q^T * rhs
} /* end namespace LeastSquares */

template <typename Scalar>
class GenLeastSquaresSolver {
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

  // Algebraic operations:
  // Complete cycle is: READY -> FACTORED -> PROJECTED -> SOLVED
  // Several steps can be accomplished in one call to statusIs()
  // Status must be READY before writing in the matrix buffer
  // Status must be READY or FACTORED before writing in the rhs buffer
  typedef LeastSquares::Status Status;
  Status status() const { return status_; }
  void statusIs(Status);

  // Ctor 
  GenLeastSquaresSolver();

private:
  static void assertOverdeterminedSystem(int eqnCount, int unknownCount);

  int equationCount_, unknownCount_;
  int largestDimension_, rhsCount_;
  int matrixLeadDim_, rhsLeadDim_;

  Status status_;

  typedef SimpleBuffer<Scalar> ScalarBuffer;
  
  ScalarBuffer matrixBuffer_;
  ScalarBuffer rhsBuffer_;
  ScalarBuffer tauBuffer_;

  // Implementation
  void factor();
  void project();
  void solve();
  void unsolve();

  // Disallow copy and assignment
  GenLeastSquaresSolver(const GenLeastSquaresSolver &);
  GenLeastSquaresSolver &operator=(const GenLeastSquaresSolver &);
};

// Inline member functions

template <typename Scalar>
inline
const Scalar *
GenLeastSquaresSolver<Scalar>::matrixBuffer() const {
  return matrixBuffer_.array();
} 

template <typename Scalar>
inline
const Scalar *
GenLeastSquaresSolver<Scalar>::matrixColBuffer(int col) const {
  return matrixBuffer() + (col * matrixLeadDim_);
}

template <typename Scalar>
inline
Scalar
GenLeastSquaresSolver<Scalar>::matrixEntry(int row, int col) const {
  return matrixColBuffer(col)[row];
}

template <typename Scalar>
inline
Scalar *
GenLeastSquaresSolver<Scalar>::matrixBuffer() {
  return matrixBuffer_.array();
} 

template <typename Scalar>
inline
Scalar *
GenLeastSquaresSolver<Scalar>::matrixColBuffer(int col) {
  return matrixBuffer() + (col * matrixLeadDim_);
}

template <typename Scalar>
inline
Scalar &
GenLeastSquaresSolver<Scalar>::matrixEntry(int row, int col) {
  return matrixColBuffer(col)[row];
}

template <typename Scalar>
inline
const Scalar *
GenLeastSquaresSolver<Scalar>::rhsBuffer() const {
  return rhsBuffer_.array();
} 

template <typename Scalar>
inline
const Scalar *
GenLeastSquaresSolver<Scalar>::rhsBuffer(int rank) const {
  return rhsBuffer() + (rank * rhsLeadDim_);
}

template <typename Scalar>
inline
Scalar
GenLeastSquaresSolver<Scalar>::rhsEntry(int rank, int row) const {
  return rhsBuffer(rank)[row];
}

template <typename Scalar>
inline
Scalar
GenLeastSquaresSolver<Scalar>::rhsEntry(int row) const {
  return rhsEntry(0, row);
}

template <typename Scalar>
inline
Scalar *
GenLeastSquaresSolver<Scalar>::rhsBuffer() {
  return rhsBuffer_.array();
} 

template <typename Scalar>
inline
Scalar *
GenLeastSquaresSolver<Scalar>::rhsBuffer(int rank) {
  return rhsBuffer() + (rank * rhsLeadDim_);
}

template <typename Scalar>
inline
Scalar &
GenLeastSquaresSolver<Scalar>::rhsEntry(int rank, int row) {
  return rhsBuffer(rank)[row];
}

template <typename Scalar>
inline
Scalar &
GenLeastSquaresSolver<Scalar>::rhsEntry(int row) {
  return rhsEntry(0, row);
}

// Implementation

template <typename Scalar>
GenLeastSquaresSolver<Scalar>::GenLeastSquaresSolver() :
  equationCount_(0),
  unknownCount_(0),
  largestDimension_(0),
  rhsCount_(0),
  matrixLeadDim_(1),
  rhsLeadDim_(1),
  status_(LeastSquares::READY)
{}

template <typename Scalar>
void
GenLeastSquaresSolver<Scalar>::problemSizeIs(int eqnCount, int unknownCount, int rhsCount) {
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

  status_ = LeastSquares::READY;
}

template <typename Scalar>
void
GenLeastSquaresSolver<Scalar>::assertOverdeterminedSystem(int eqnCount, int unknownCount) {
  if (eqnCount < unknownCount) {
    throw std::domain_error("Underdetermined system");
  }
}

} /* end namespace Rom */

#endif /* ROM_LEASTSQUARESSOLVER_H */
