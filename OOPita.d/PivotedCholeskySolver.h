#ifndef PITA_PIVOTEDCHOLESKYSOLVER_H
#define PITA_PIVOTEDCHOLESKYSOLVER_H

#include "Fwk.h"

#include <Math.d/Vector.h>
#include <Math.d/FullSquareMatrix.h>
#include <Math.d/SymFullMatrix.h>
#include "SimpleBuffer.h"

namespace Pita {

class PivotedCholeskySolver : public Fwk::PtrInterface<PivotedCholeskySolver> {
public:
  typedef Fwk::Ptr<PivotedCholeskySolver> Ptr;
  typedef Fwk::Ptr<const PivotedCholeskySolver> PtrConst;

  enum Status {
    NON_FACTORIZED,
    FACTORIZED
  };

  Status status() const { return status_; }
  int matrixSize() const { return matrixSize_; }

  // These 3 accessors have meaningful values only if status() == FACTORIZED
  int factorRank() const { return factorRank_; }
  int factorPermutation(int index) const { return factorPermutation_[index] - 1; }
  const FullSquareMatrix & choleskyFactor() const { return choleskyFactor_; }

  // Mutators
  void matrixIs(const SymFullMatrix & matrix);
  void statusIs(Status s);
  
  // Solve in place: Input parameter rhs is modified
  Vector & solution(Vector & rhs) const;

  static Ptr New() {
    return new PivotedCholeskySolver();
  }

protected:
  PivotedCholeskySolver();

private:
  Status status_;
  int matrixSize_;
  int factorRank_;
  FullSquareMatrix choleskyFactor_;
  SimpleBuffer<int> factorPermutation_; // Indices starting at 1 (Fortran convention)

  DISALLOW_COPY_AND_ASSIGN(PivotedCholeskySolver);
}; 

} // end namespace Pita

#endif /* PITA_PIVOTEDCHOLESKYSOLVER_H */
