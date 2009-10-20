#ifndef PITA_PIVOTEDCHOLESKYSOLVER_H
#define PITA_PIVOTEDCHOLESKYSOLVER_H

#include "Fwk.h"
#include "RankDeficientSolver.h"
#include <Math.d/Vector.h>
#include <Math.d/FullSquareMatrix.h>
#include <Math.d/SymFullMatrix.h>
#include "SimpleBuffer.h"

namespace Pita {

class PivotedCholeskySolver : public RankDeficientSolver {
public:
  EXPORT_PTRINTERFACE_TYPES(PivotedCholeskySolver);

  // Accessors
  const FullSquareMatrix & choleskyFactor() const { return choleskyFactor_; }
  virtual const Vector & solution(Vector & rhs) const; // In-place solution: rhs modified

  // Mutators
  void matrixIs(const SymFullMatrix & matrix);
  void matrixIs(const FullSquareMatrix & matrix); // Use only lower triangular part
  virtual void transposedMatrixIs(const FullSquareMatrix & matrix) { matrixIs(matrix); } // Use only upper triangular part
  virtual void orderingIs(Ordering o);
  
  static Ptr New(double tolerance = -1.0) {
    return new PivotedCholeskySolver(tolerance);
  }

protected:
  explicit PivotedCholeskySolver(double tolerance);

  void performFactorization();

private:
  FullSquareMatrix choleskyFactor_;

  DISALLOW_COPY_AND_ASSIGN(PivotedCholeskySolver);
}; 

} // end namespace Pita

#endif /* PITA_PIVOTEDCHOLESKYSOLVER_H */
