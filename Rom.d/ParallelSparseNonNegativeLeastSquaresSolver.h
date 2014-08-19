#ifndef ROM_PARALLELSPARSENONNEGATIVELEASTSQUARESSOLVER_H
#define ROM_PARALLELSPARSENONNEGATIVELEASTSQUARESSOLVER_H

#include "SparseNonNegativeLeastSquaresSolver.h"
#include "SimpleBuffer.h"
#include <vector>

namespace Rom {

class ParallelSparseNonNegativeLeastSquaresSolver {
public:

  int subdomainCount() const { return nsub_; }
  SparseNonNegativeLeastSquaresSolver<std::vector<double>,size_t> * subdomainSolver(int i) { return sd_[i]; }

  // Problem size
  long equationCount() const { return equationCount_; }
  long unknownCount() const { return unknownCount_; }

  void problemSizeIs(long eqnCount, long unkCount);

  // Options
  double relativeTolerance() const { return relativeTolerance_; }
  void relativeToleranceIs(double relTol) { relativeTolerance_ = relTol; }

  bool verboseFlag() const { return verboseFlag_; }
  void verboseFlagIs(bool verFlg) { verboseFlag_ = verFlg; }

  bool scalingFlag() const { return scalingFlag_; }
  void scalingFlagIs(bool scaFlg) { scalingFlag_ = scaFlg; }

  bool positivityFlag() const { return positivity_; }
  void positivityIs(bool posFlg) { positivity_ = posFlg; }

  int solverType() const { return solverType_; }
  void solverTypeIs(int solTyp) { solverType_ = solTyp; }

  double maxSizeRatio() const { return maxSizeRatio_; }
  void maxSizeRatioIs(double maxSze) { maxSizeRatio_ = maxSze; }

  double maxIterRatio() const { return maxIterRatio_; }
  void maxIterRatioIs(double maxIte) { maxIterRatio_ = maxIte; }

  // Rhs buffer: [equationCount]
  const double * rhsBuffer() const { return rhsBuffer_.array(); }
  double * rhsBuffer() { return rhsBuffer_.array(); }

  // Error magnitude
  double errorMagnitude() const { return errorMagnitude_; }

  // Perform solve
  void solve();

  // Constructor
  ParallelSparseNonNegativeLeastSquaresSolver(int nsub, SparseNonNegativeLeastSquaresSolver<std::vector<double>,size_t> **sd);

private:

  long equationCount_;
  long unknownCount_;

  SimpleBuffer<double> rhsBuffer_;

  double relativeTolerance_;
  double errorMagnitude_;
  bool verboseFlag_;
  bool scalingFlag_;
  bool positivity_;
  int solverType_;
  double maxSizeRatio_;
  double maxIterRatio_;

  int nsub_;
  SparseNonNegativeLeastSquaresSolver<std::vector<double>,size_t> **sd_;

  // Disallow copy & assignment
  ParallelSparseNonNegativeLeastSquaresSolver(const ParallelSparseNonNegativeLeastSquaresSolver &);
  ParallelSparseNonNegativeLeastSquaresSolver &operator=(const ParallelSparseNonNegativeLeastSquaresSolver &);
};

} // end namespace Rom

#endif /* ROM_PARALLELSPARSENONNEGATIVELEASTSQUARESSOLVER_H */

