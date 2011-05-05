#ifndef ROM_GAUSSNEWTONSOLVER_H
#define ROM_GAUSSNEWTONSOLVER_H

#include "PodProjectionSolver.h"
#include "LeastSquaresSolver.h"

#include <cassert>
#include <algorithm>

namespace Rom {

template <typename Scalar>
class GenGaussNewtonSolver : public GenPodProjectionSolver<Scalar> {
public:
  GenGaussNewtonSolver(Connectivity *cn, DofSetArray *dsa, ConstrainedDSA *c_dsa);
  
private:
  GenLeastSquaresSolver<Scalar> lsSolver_;

  // Overriden implementation functions
  virtual void resetSolver(int vCount, int vSize);
  virtual void assembleAndFactorReducedSystem();
  virtual void projectRhs(const GenVector<Scalar> &);
  virtual double getReducedRhsNorm() const;
  virtual void solveReducedSystem(GenVector<Scalar> &);

  void fillRhsBuffer(const GenVector<Scalar> &);
};

template <typename Scalar>
GenGaussNewtonSolver<Scalar>::GenGaussNewtonSolver(Connectivity *cn,
                                                   DofSetArray *dsa,
                                                   ConstrainedDSA *c_dsa):
  GenPodProjectionSolver<Scalar>(cn, dsa, c_dsa),
  lsSolver_()
{}

template <typename Scalar>
void
GenGaussNewtonSolver<Scalar>::assembleAndFactorReducedSystem() {
  lsSolver_.statusIs(LeastSquares::READY);
  
  const Scalar* const matrixBufferBegin = &(this->lastReducedMatrixAction()[0][0]);
  std::copy(matrixBufferBegin,
            matrixBufferBegin + lsSolver_.equationCount() * lsSolver_.unknownCount(),
            lsSolver_.matrixBuffer());

  lsSolver_.statusIs(LeastSquares::FACTORED);
}

template <typename Scalar>
inline
void
GenGaussNewtonSolver<Scalar>::fillRhsBuffer(const GenVector<Scalar> &rhs) {
  std::copy(rhs.data(), rhs.data() + rhs.size(), lsSolver_.rhsBuffer(0));
}

template <typename Scalar>
void
GenGaussNewtonSolver<Scalar>::resetSolver(int vCount, int vSize) {
  lsSolver_.problemSizeIs(vSize, vCount);
}

template <typename Scalar>
void
GenGaussNewtonSolver<Scalar>::projectRhs(const GenVector<Scalar> &rhs) {
  assert(lsSolver_.status() == LeastSquares::FACTORED);

  fillRhsBuffer(rhs);
  lsSolver_.statusIs(LeastSquares::PROJECTED);
}

template <typename Scalar>
void
GenGaussNewtonSolver<Scalar>::solveReducedSystem(GenVector<Scalar> &rhs) {
  assert(lsSolver_.status() == LeastSquares::FACTORED ||
         lsSolver_.status() == LeastSquares::PROJECTED);
  
  if (lsSolver_.status() == LeastSquares::FACTORED) {
    fillRhsBuffer(rhs);
  }

  lsSolver_.statusIs(LeastSquares::SOLVED);
  std::copy(lsSolver_.rhsBuffer(0),
            lsSolver_.rhsBuffer(0) + lsSolver_.unknownCount(),
            this->getReducedSolution().data());
  lsSolver_.statusIs(LeastSquares::FACTORED);
}

template <typename Scalar>
double
GenGaussNewtonSolver<Scalar>::getReducedRhsNorm() const {
  const GenVector<Scalar> components(const_cast<Scalar *>(lsSolver_.rhsBuffer(0)),
                                     lsSolver_.unknownCount(),
                                     false);
  return components.norm();
}

typedef GenGaussNewtonSolver<double> GaussNewtonSolver;

} /* end namespace Rom */

#endif /* ROM_GAUSSNEWTONSOLVER_H */
