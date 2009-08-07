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

  const FullSquareMatrix & choleskyFactor() const { return choleskyFactor_; }
  double tolerance() const { return tolerance_; } // -1.0 stands for default tolerance

  // Mutators
  void matrixIs(const SymFullMatrix & matrix);
  void matrixIs(const FullSquareMatrix & matrix); // Use only lower triangular part
  void statusIs(Status s);
  
  virtual const Vector & solution(Vector & rhs) const; // In-place solution: rhs modified

  static Ptr New(double tolerance = -1.0) {
    return new PivotedCholeskySolver(tolerance);
  }

protected:
  explicit PivotedCholeskySolver(double tolerance);

private:
  FullSquareMatrix choleskyFactor_;
  double tolerance_;

  DISALLOW_COPY_AND_ASSIGN(PivotedCholeskySolver);
}; 

} // end namespace Pita

#endif /* PITA_PIVOTEDCHOLESKYSOLVER_H */
