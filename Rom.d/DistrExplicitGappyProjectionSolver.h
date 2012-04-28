#ifndef ROM_DISTREXPLICITGAPPYPROJECTIONSOLVER_H
#define ROM_DISTREXPLICITGAPPYPROJECTIONSOLVER_H

#include <Solvers.d/ParallelSolver.h>

#include "DistrVecBasis.h"
#include "DistrVecBasisOps.h"

#include <Timers.d/Timing.h>
#include <Paral.d/SubDOp.h>

#include <Math.d/Vector.h>
#include <Math.d/FullSquareMatrix.h>

#include <algorithm>
#include <memory>
#include <cassert>

namespace Rom {

template <typename Scalar>
class GenDistrExplicitGappyProjectionSolver : public GenParallelSolver<Scalar> {
public:
  virtual void refactor();
  virtual void solve(GenDistrVector<Scalar> &, GenDistrVector<Scalar> &);
  virtual void reSolve(GenDistrVector<Scalar> &);
 
  virtual int numRBM() { return 0; }
  virtual void getRBMs(Scalar *) {}
  virtual void getRBMs(GenDistrVectorSet<Scalar> &) {}
  
  virtual Timings &getTimers() { return timers_; }
  virtual double getSolutionTime() { return 0.0; }

  // Reduced basis parameters (read-only)
  typedef GenVecBasis<Scalar, GenDistrVector> BasisType;
  const BasisType &reconstructionBasis() const { return *reconstructionBasis_; }
  const BasisType &projectionBasis() const { return *projectionBasis_; }
  
  // Gappy operators
  // Passed bases must be kept alive by owner
  explicit GenDistrExplicitGappyProjectionSolver(const BasisType &projBasis, const BasisType &reconstrBasis);

  virtual ~GenDistrExplicitGappyProjectionSolver();
 
private:
  Timings timers_;

  const BasisType *reconstructionBasis_;
  const BasisType *projectionBasis_;
  
  GenVector<Scalar> components_; 

  // Disallow copy & assignment
  GenDistrExplicitGappyProjectionSolver(const GenDistrExplicitGappyProjectionSolver &);
  GenDistrExplicitGappyProjectionSolver &operator=(const GenDistrExplicitGappyProjectionSolver &);
};

template <typename Scalar>
GenDistrExplicitGappyProjectionSolver<Scalar>::GenDistrExplicitGappyProjectionSolver(const BasisType &projBasis, const BasisType &reconstrBasis) :
  timers_(),
  reconstructionBasis_(&reconstrBasis),
  projectionBasis_(&projBasis),
  components_(projBasis.vectorCount())
{
  assert(projectionBasis_->vectorCount() == reconstructionBasis_->vectorCount());
}

template <typename Scalar>
GenDistrExplicitGappyProjectionSolver<Scalar>::~GenDistrExplicitGappyProjectionSolver() {
  // Nothing to do
}

template <typename Scalar>
void
GenDistrExplicitGappyProjectionSolver<Scalar>::refactor() {
  // Nothing to do (Solver is entirely passive)
}

template <typename Scalar>
void
GenDistrExplicitGappyProjectionSolver<Scalar>::solve(GenDistrVector<Scalar> &rhs, GenDistrVector<Scalar> &result) {
  // TODO: Restricted projection
  vector_components(*projectionBasis_, rhs, components_);
  assembled_vector(*reconstructionBasis_, components_, result);
}

template <typename Scalar>
void
GenDistrExplicitGappyProjectionSolver<Scalar>::reSolve(GenDistrVector<Scalar> &rhs) {
  solve(rhs, rhs);
}

typedef GenDistrExplicitGappyProjectionSolver<double> DistrExplicitGappyProjectionSolver;

} // end namespace Rom

#endif /* ROM_DISTREXPLICITGAPPYPROJECTIONSOLVER_H */
