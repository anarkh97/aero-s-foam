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
  virtual void refactor();
  virtual void solve(GenDistrVector<Scalar> &, GenDistrVector<Scalar> &);
  virtual void reSolve(GenDistrVector<Scalar> &);
 
  virtual int numRBM() { return 0; }
  virtual void getRBMs(Scalar *) {}
  virtual void getRBMs(GenDistrVectorSet<Scalar> &) {}
  
  virtual Timings &getTimers() { return timers_; }
  virtual double getSolutionTime() { return 0.0; }

  // Full-order matrix (assumed symmetric positive definite)
  const GenSubDOp<Scalar> &fullMatrix() const { return fullMatrix; }

  // Reduced basis parameters
  typedef GenVecBasis<Scalar, GenDistrVector> BasisType;
  const BasisType &projectionBasis() const { return *projectionBasis_; }
  void projectionBasisIs(const BasisType &); // Passed object must be kept alive by owner

  explicit GenDistrGalerkinProjectionSolver(const GenSubDOp<Scalar> &fullMat); // Passed object must be kept alive by owner

  virtual ~GenDistrGalerkinProjectionSolver();

  double counter;
  double timer;

private:
  Timings timers_;

  const GenSubDOp<Scalar> *fullMatrix_;

  const BasisType *projectionBasis_; // Assumed to hold fully consistent vectors
  BasisType normalizedBasis_;

  // Disallow copy & assignment
  GenDistrGalerkinProjectionSolver(const GenDistrGalerkinProjectionSolver &);
  GenDistrGalerkinProjectionSolver &operator=(const GenDistrGalerkinProjectionSolver &);
};

template <typename Scalar>
GenDistrGalerkinProjectionSolver<Scalar>::GenDistrGalerkinProjectionSolver(const GenSubDOp<Scalar> &fullMat) :
  timers_(),
  fullMatrix_(&fullMat),
  normalizedBasis_()
{
  projectionBasis_ = &normalizedBasis_;
}

template <typename Scalar>
GenDistrGalerkinProjectionSolver<Scalar>::~GenDistrGalerkinProjectionSolver() {
  // Nothing to do
}

template <typename Scalar>
void
GenDistrGalerkinProjectionSolver<Scalar>::projectionBasisIs(const BasisType &reducedBasis) {
  projectionBasis_ = &reducedBasis;
  refactor();
}

template <typename Scalar>
void
GenDistrGalerkinProjectionSolver<Scalar>::refactor() {
  renormalized_basis(*fullMatrix_, *projectionBasis_, normalizedBasis_);
  timer = 0;
  counter = 0;
}

template <typename Scalar>
void
GenDistrGalerkinProjectionSolver<Scalar>::solve(GenDistrVector<Scalar> &rhs, GenDistrVector<Scalar> &result) {
  const int vectorCount = normalizedBasis_.vectorCount();
  GenVector<Scalar> components(vectorCount, Scalar()); 

  counter += 1;
  double dummy = 0;
  
  if(domain->solInfo().elemLumpPodRom){

    dummy -= getTime();
    normalizedBasis_.project(rhs, result);  
    dummy += getTime();
    timer += dummy;
    std::cout << "             modelIII projection time = " << timer/counter << std::endl;

  } else {

    dummy -= getTime();
    //---------------------------------------------------------------------
    vector_components_vector_masterflag(normalizedBasis_, rhs, components);
    //---------------------------------------------------------------------
    assembled_vector(normalizedBasis_, components, result);
    //-----------------------------------------------------
    dummy += getTime();
    timer += dummy;
    std::cout << "             modelII  projection time = " << timer/counter << std::endl;

       }
}

template <typename Scalar>
void
GenDistrGalerkinProjectionSolver<Scalar>::reSolve(GenDistrVector<Scalar> &rhs) {
  solve(rhs, rhs);
}

typedef GenDistrGalerkinProjectionSolver<double> DistrGalerkinProjectionSolver;

} // end namespace Rom

#endif /* ROM_DISTRGALERKINPROJECTIONSOLVER_H */
