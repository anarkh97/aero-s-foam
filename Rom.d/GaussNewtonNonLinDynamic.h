#ifndef ROM_GAUSSNEWTONNONLINDYNAMIC_H
#define ROM_GAUSSNEWTONNONLINDYNAMIC_H

#include "GalerkinProjectionSolver.h"

#include <Problems.d/NonLinDynam.h>
#include <Math.d/VectorSet.h>

#include <memory>

class GaussNewtonNonLinDynamic : public NonLinDynamic {
public:
  explicit GaussNewtonNonLinDynamic(Domain *);

  // Required additional pre- and post-processing
  virtual void preProcess();
  void postProcess();

  GalerkinProjectionSolver *getSolver(); // Hiding NonLinDynamic::getSolver

  // Helper class to be used as template parameter in NLDynamSolver 
  class Updater;
 
private:
  virtual bool factorWhenBuilding() const; // Overriden

  std::auto_ptr<GenVectorSet<double> > projectionBasis_;
};

#endif /* ROM_GAUSSNEWTONNONLINDYNAMIC_H */
