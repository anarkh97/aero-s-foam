#ifndef ROM_PODPROJECTIONNONLINDYNAMIC_H
#define ROM_PODPROJECTIONNONLINDYNAMIC_H

#include "PodProjectionSolver.h"
#include "VecNodeDof6Conversion.h"
#include "NodeDof6Buffer.h"
#include "BasisFileStream.h"

#include <Problems.d/NonLinDynam.h>
#include <Rom.d/VecBasis.h>

#include <memory>

namespace Rom {

class PodProjectionNonLinDynamic : public NonLinDynamic {
public:
  explicit PodProjectionNonLinDynamic(Domain *);

  // Required additional pre-processing
  virtual void preProcess();

  // Hiding NonLinDynamic::getSolve
  PodProjectionSolver *getSolver();
  const PodProjectionSolver *getSolver() const;

  // Helper class to be used as template parameter in NLDynamSolver 
  class Updater;

  // Hooks in NLDynamSolver
  virtual double getResidualNorm(const Vector &);
  int checkConvergence(int iteration, double normRes, Vector &residual, Vector &dv, double time); // relies on function hiding

private:
  virtual bool factorWhenBuilding() const; // Overriden

  void computeAndSaveJacobianSnapshot();
  void saveResidualSnapshot(const Vector &snap) { *residualSnapFile_ << snap; }

  VecBasis projectionBasis_;

  std::auto_ptr<VecNodeDof6Conversion> vecNodeDof6Conversion_;
  NodeDof6Buffer snapBuffer_;
  
  std::auto_ptr<BasisOutputStream> residualSnapFile_;
  std::auto_ptr<BasisOutputStream> jacobianSnapFile_;

  friend class Updater;

  // Disallow copy and assignment
  PodProjectionNonLinDynamic(const PodProjectionNonLinDynamic &);
  PodProjectionNonLinDynamic &operator=(const PodProjectionNonLinDynamic &);
};

#include <Driver.d/StateUpdater.h>

// Provides hooks to be used in NLDynamSolver to call the snapshot collection functions
class PodProjectionNonLinDynamic::Updater : public IncrUpdater<PodProjectionNonLinDynamic, GenVector<double>, GeomState> {
public:
  static double formRHScorrector(PodProjectionNonLinDynamic *pbd, GenVector<double> &inc_displac,
                                 GenVector<double> &vel_n, GenVector<double> &accel,
                                 GenVector<double> &residual, GenVector<double> &rhs,
                                 GeomState *geomState) {
    const double result = IncrUpdater<PodProjectionNonLinDynamic, GenVector<double>, GeomState>::formRHScorrector(
        pbd, inc_displac, vel_n, accel, residual, rhs, geomState);

    pbd->saveResidualSnapshot(rhs);

    return result;
  }
};

} /* end namespace Rom */

#endif /* ROM_PODPROJECTIONNONLINDYNAMIC_H */
