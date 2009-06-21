#ifndef PITA_NEARSYMMETRICSOLVER_H
#define PITA_NEARSYMMETRICSOLVER_H

#include "Fwk.h"
#include "RankDeficientSolver.h"
#include <Math.d/FullSquareMatrix.h>

namespace Pita {

class NearSymmetricSolver : public RankDeficientSolver {
public:
  EXPORT_PTRINTERFACE_TYPES(NearSymmetricSolver);

  double tolerance() const { return tolerance_; }
  const FullSquareMatrix & transposedMatrix() const { return transposedMatrix_; }

  void toleranceIs(double tol);
  void transposedMatrixIs(const FullSquareMatrix & tm); 
  void statusIs(Status s);

  virtual const Vector & solution(Vector & rhs) const; // In-place solution: rhs modified

  static Ptr New(double tolerance) {
    return new NearSymmetricSolver(tolerance);
  }

protected:
  explicit NearSymmetricSolver(double tol);

  void setTolerance(double tol) { tolerance_ = tol; }
  void setTransposedMatrix(const FullSquareMatrix & tm) { transposedMatrix_ = tm; } 

private:
  double tolerance_;
  FullSquareMatrix transposedMatrix_;

  SimpleBuffer<int> pivots_;
};

} /* end namespace Pita */

#endif /* PITA_NEARSYMMETRICSOLVER_H */
