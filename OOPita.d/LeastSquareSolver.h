#ifndef PITA_LEASTSQUARESOLVER_H
#define PITA_LEASTSQUARESOLVER_H

#include "Fwk.h"

#include <Math.d/FullSquareMatrix.h>
#include <Math.d/Vector.h>

#include "SimpleBuffer.h"

namespace Pita {

class LeastSquareSolver : public Fwk::PtrInterface<LeastSquareSolver> {
public:
  EXPORT_PTRINTERFACE_TYPES(LeastSquareSolver);
  
  enum Status {
    NON_FACTORIZED,
    FACTORIZED
  };

  Status status() const { return status_; }
  double tolerance() const { return tolerance_; }
  int matrixSize() const { return matrixSize_; }
  const FullSquareMatrix & transposedMatrix() const { return transposedMatrix_; }
  
  // These 2 accessors have meaningful values only if status() == FACTORIZED
  int factorRank() const { return factorRank_; }
  int factorPermutation(int index) const { return factorPermutation_[index] - 1; }

  void transposedMatrixIs(const FullSquareMatrix & transposedMatrix); // The matrix is copied
  void statusIs(Status s);
  void toleranceIs(double tol);

  // In place solution, argument rhs is modified
  const Vector & solution(Vector & rhs) const;

  static Ptr New(double tolerance) {
    return new LeastSquareSolver(tolerance);
  }

protected:
  explicit LeastSquareSolver(double tol);
  
  void updateFactorRank();

private:
  Status status_;
  FullSquareMatrix transposedMatrix_; // To have column-major ordering
  int matrixSize_;
  double tolerance_;
  int factorRank_;

  SimpleBuffer<int> factorPermutation_;
  SimpleBuffer<double> tau_;
  mutable SimpleBuffer<double> workspace_;

  DISALLOW_COPY_AND_ASSIGN(LeastSquareSolver);
};

} /* end namespace Pita */

#endif /* PITA_LEASTSQUARESOLVER_H */
