#ifndef ROM_GAUSSNEWTONNONLINDYNAMIC_H
#define ROM_GAUSSNEWTONNONLINDYNAMIC_H

#include "GaussNewtonSolver.h"
#include "VecNodeDof6Conversion.h"
#include "NodeDof6Buffer.h"
#include "BasisFileStream.h"

#include <Problems.d/NonLinDynam.h>
#include <Rom.d/VecBasis.h>

#include <memory>

class GaussNewtonNonLinDynamic : public NonLinDynamic {
public:
  explicit GaussNewtonNonLinDynamic(Domain *);

  // Required additional pre-processing
  virtual void preProcess();

  // Hiding NonLinDynamic::getSolve
  GaussNewtonSolver *getSolver();
  const GaussNewtonSolver *getSolver() const;

  // Helper class to be used as template parameter in NLDynamSolver 
  class Updater;

  // Direct hooks in NLDynamSolver (rely on function hiding)
  int checkConvergence(int iteration, double normRes, Vector &residual, Vector &dv, double time);
  double getResidualNorm(const Vector &);

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
  GaussNewtonNonLinDynamic(const GaussNewtonNonLinDynamic &);
  GaussNewtonNonLinDynamic &operator=(const GaussNewtonNonLinDynamic &);
};

#include <Driver.d/StateUpdater.h>

// Provides hooks to be used in NLDynamSolver to call the snapshot collection functions
class GaussNewtonNonLinDynamic::Updater : public IncrUpdater<GaussNewtonNonLinDynamic, GenVector<double>, GeomState> {
public:
  static double formRHScorrector(GaussNewtonNonLinDynamic *pbd, GenVector<double> &inc_displac,
                                 GenVector<double> &vel_n, GenVector<double> &accel,
                                 GenVector<double> &residual, GenVector<double> &rhs,
                                 GeomState *geomState) {
    const double result = IncrUpdater<GaussNewtonNonLinDynamic, GenVector<double>, GeomState>::formRHScorrector(
        pbd, inc_displac, vel_n, accel, residual, rhs, geomState);

    pbd->saveResidualSnapshot(rhs);

    return result;
  }
};
#endif /* ROM_GAUSSNEWTONNONLINDYNAMIC_H */
