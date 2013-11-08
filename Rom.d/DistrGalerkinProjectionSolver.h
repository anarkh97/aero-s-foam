#ifndef ROM_DISTRGALERKINPROJECTIONSOLVER_H
#define ROM_DISTRGALERKINPROJECTIONSOLVER_H

#include <Solvers.d/ParallelSolver.h>

#include "DistrVecBasis.h"
#include "DistrVecBasisOps.h"

#include <Timers.d/Timing.h>
#include <Paral.d/SubDOp.h>

#include <Math.d/Vector.h>
#include <Math.d/FullSquareMatrix.h>

#include <algorithm>
#include <memory>

namespace Rom {

template <typename Scalar>
class GenDistrGalerkinProjectionSolver : public GenParallelSolver<Scalar> {
public:
  virtual void solve(GenDistrVector<Scalar> &, GenDistrVector<Scalar> &);
  virtual void reSolve(GenDistrVector<Scalar> &);
 
  virtual int numRBM() { return 0; }
  virtual void getRBMs(Scalar *) {}
  virtual void getRBMs(GenDistrVectorSet<Scalar> &) {}
  
  virtual Timings &getTimers() { return timers_; }
  virtual double getSolutionTime() { return 0.0; }

  // Reduced basis parameters
  typedef GenVecBasis<Scalar, GenDistrVector> BasisType;

  explicit GenDistrGalerkinProjectionSolver(const BasisType &); // Passed object must be kept alive by owner

  virtual ~GenDistrGalerkinProjectionSolver();

  double counter;
  double timer;

private:
  Timings timers_;

  BasisType normalizedBasis_;

  // Disallow copy & assignment
  GenDistrGalerkinProjectionSolver(const GenDistrGalerkinProjectionSolver &);
  GenDistrGalerkinProjectionSolver &operator=(const GenDistrGalerkinProjectionSolver &);
};

template <typename Scalar>
GenDistrGalerkinProjectionSolver<Scalar>::GenDistrGalerkinProjectionSolver(const BasisType &NormBasis) :
  timers_(),
  normalizedBasis_(NormBasis)
{}

template <typename Scalar>
GenDistrGalerkinProjectionSolver<Scalar>::~GenDistrGalerkinProjectionSolver() {
  // Nothing to do
}

template <typename Scalar>
void
GenDistrGalerkinProjectionSolver<Scalar>::solve(GenDistrVector<Scalar> &rhs, GenDistrVector<Scalar> &result) {
  result = rhs;
}

template <typename Scalar>
void
GenDistrGalerkinProjectionSolver<Scalar>::reSolve(GenDistrVector<Scalar> &rhs) {
  // Nothing to do
}

typedef GenDistrGalerkinProjectionSolver<double> DistrGalerkinProjectionSolver;

} // end namespace Rom

#endif /* ROM_DISTRGALERKINPROJECTIONSOLVER_H */
