#ifndef PITA_REDUCEDFULLTIMESLICEIMPL_H
#define PITA_REDUCEDFULLTIMESLICEIMPL_H

#include "ReducedFullTimeSlice.h"

#include "RankDeficientSolver.h"

template <typename Scalar> class GenFullSquareMatrix;
typedef GenFullSquareMatrix<double> FullSquareMatrix;

namespace Pita {

class ReducedFullTimeSliceImpl : public ReducedFullTimeSlice {
public:
  EXPORT_PTRINTERFACE_TYPES(ReducedFullTimeSliceImpl);
  class Manager;

  // Overriden mutators
  virtual void iterationIs(IterationRank i);

  // Reduced-space operators
  const FullSquareMatrix * reprojectionMatrix() const { return reprojectionMatrix_; }
  const RankDeficientSolver * projectionSolver() const { return projectionSolver_; }

  void reprojectionMatrixIs(const FullSquareMatrix * rm) { reprojectionMatrix_ = rm; }
  void projectionSolverIs(const RankDeficientSolver * ps) { projectionSolver_ = ps; }

protected:
  ReducedFullTimeSliceImpl(HalfSliceRank headRank,
                           const FullSquareMatrix * reprojectionMatrix,
                           const RankDeficientSolver * projectionSolver);
  
private:
  const FullSquareMatrix * reprojectionMatrix_;
  const RankDeficientSolver * projectionSolver_;
};


class ReducedFullTimeSliceImpl::Manager : public ReducedFullTimeSlice::Manager, private Fwk::GenManagerImpl<ReducedFullTimeSliceImpl, HalfSliceRank> {
public:
  EXPORT_PTRINTERFACE_TYPES(Manager);

  ReducedFullTimeSliceImpl * instance(const HalfSliceRank & r) const { return Fwk::GenManagerImpl<ReducedFullTimeSliceImpl, HalfSliceRank>::instance(r); }
  InstanceCount instanceCount() const { return Fwk::GenManagerImpl<ReducedFullTimeSliceImpl, HalfSliceRank>::instanceCount(); }

  ReducedFullTimeSliceImpl * instanceNew(const HalfSliceRank & r) { return Fwk::GenManagerImpl<ReducedFullTimeSliceImpl, HalfSliceRank>::instanceNew(r); } 
  void instanceDel(const HalfSliceRank & r) { Fwk::GenManagerImpl<ReducedFullTimeSliceImpl, HalfSliceRank>::instanceDel(r); }

  const FullSquareMatrix * defaultReprojectionMatrix() const { return defaultReprojectionMatrix_; }
  const RankDeficientSolver * defaultProjectionSolver() const { return defaultProjectionSolver_.ptr(); }
  
  void defaultReprojectionMatrixIs(const FullSquareMatrix * drm) { defaultReprojectionMatrix_ = drm; }
  void defaultPojectionSolverIs(const RankDeficientSolver * ds) { defaultProjectionSolver_ = ds; }

  static Ptr New(const FullSquareMatrix * defaultReprojectionMatrix, const RankDeficientSolver * defaultProjectionSolver) {
    return new Manager(defaultReprojectionMatrix, defaultProjectionSolver);
  }
  
protected:
  Manager(const FullSquareMatrix * defaultReprojectionMatrix, const RankDeficientSolver * defaultProjectionSolver);

  virtual ReducedFullTimeSliceImpl * createNewInstance(const HalfSliceRank & r);

private:
  const FullSquareMatrix * defaultReprojectionMatrix_;
  RankDeficientSolver::PtrConst defaultProjectionSolver_;
};

} // end namespace Pita

#endif /* PITA_REDUCEDFULLTIMESLICEIMPL_H */
