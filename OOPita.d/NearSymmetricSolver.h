#ifndef PITA_NEARSYMMETRICSOLVER_H
#define PITA_NEARSYMMETRICSOLVER_H

#include "Fwk.h"
#include "RankDeficientSolver.h"
#include <Math.d/FullSquareMatrix.h>

namespace Pita {

class NearSymmetricSolver : public RankDeficientSolver {
public:
  EXPORT_PTRINTERFACE_TYPES(NearSymmetricSolver);

  enum RescalingStatus {
    NO_RESCALING,
    ROW_RESCALING,
    SYMMETRIC_RESCALING
  };

  const FullSquareMatrix & transposedMatrix() const { return transposedMatrix_; }

  virtual void toleranceIs(double tol);
  virtual void transposedMatrixIs(const FullSquareMatrix & tm); 
  virtual void statusIs(Status s);

  virtual const Vector & solution(Vector & rhs) const; // In-place solution: rhs modified

  static Ptr New(double tolerance) {
    return new NearSymmetricSolver(tolerance);
  }

protected:
  explicit NearSymmetricSolver(double tol);

  void setTransposedMatrix(const FullSquareMatrix & tm) { transposedMatrix_ = tm; } 

private:
  FullSquareMatrix transposedMatrix_;

  SimpleBuffer<int> pivots_;

  RescalingStatus rescalingStatus_;
  SimpleBuffer<double> scaling_;
};

} /* end namespace Pita */

#endif /* PITA_NEARSYMMETRICSOLVER_H */
