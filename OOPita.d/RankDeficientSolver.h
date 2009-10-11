#ifndef PITA_RANKDEFICIENTSOLVER_H
#define PITA_RANKDEFICIENTSOLVER_H

#include "Fwk.h"
#include "SimpleBuffer.h"
#include <Math.d/Vector.h>
#include <Math.d/FullSquareMatrix.h>

namespace Pita {

class RankDeficientSolver : public Fwk::PtrInterface<RankDeficientSolver> {
public:
  EXPORT_PTRINTERFACE_TYPES(RankDeficientSolver);

  enum Status {
    NON_FACTORIZED,
    FACTORIZED
  };

  Status status() const { return status_; }
  int matrixSize() const { return matrixSize_; }
  double tolerance() const { return tolerance_; }
  int factorRank() const { return factorRank_; }
  int factorPermutation(int index) const { return factorPermutation_[index] - 1; } // 0 <= index <= factorRank

  virtual const Vector & solution(Vector & rhs) const = 0; // In-place solution: rhs modified

  virtual void toleranceIs(double tol) { setTolerance(tol); }
  virtual void transposedMatrixIs(const FullSquareMatrix & tm) = 0;
  virtual void statusIs(Status s) = 0;

protected:
  RankDeficientSolver(double tol) :
    status_(NON_FACTORIZED),
    matrixSize_(0),
    tolerance_(tol),
    factorRank_(0),
    factorPermutation_()
  {}

  const int & getFactorRank() const { return factorRank_; }
  const int & getMatrixSize() const { return matrixSize_; }
  const double & getTolerance() const { return tolerance_; }
  const SimpleBuffer<int> & getFactorPermutation() const { return factorPermutation_; }
  
  SimpleBuffer<int> & getFactorPermutation() { return factorPermutation_; }
  void setStatus(Status s) { status_ = s; }
  void setMatrixSize(int ms) { matrixSize_ = ms; }
  void setTolerance(double tol) { tolerance_ = tol; }
  void setFactorRank(int fr) { factorRank_ = fr; }

private:
  Status status_;
  int matrixSize_;
  double tolerance_;
  int factorRank_;
  SimpleBuffer<int> factorPermutation_; // Indices starting at 1 (Fortran convention)
  
  DISALLOW_COPY_AND_ASSIGN(RankDeficientSolver);
};

} /* end namespace Pita */

#endif /* PITA_RANKDEFICIENTSOLVER_H */
