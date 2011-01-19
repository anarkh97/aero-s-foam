#include "GaussNewtonNonLinDynamic.h"

#include <stdexcept>

GaussNewtonNonLinDynamic::GaussNewtonNonLinDynamic(Domain *d) :
  NonLinDynamic(d)
{}

void
GaussNewtonNonLinDynamic::preProcess() {
  NonLinDynamic::preProcess();
  
  if (!dynamic_cast<GalerkinProjectionSolver *>(NonLinDynamic::getSolver())) {
    throw std::runtime_error("Solver must be a GalerkinProjectionSolver");
  }

  // Testing hack
  projectionBasis_.reset(new GenVectorSet<double>(solVecInfo(), solVecInfo(), 0.0));
  for (int i = 0; i < solVecInfo(); ++i) {
    (*projectionBasis_)[i][i] = 1.0;
  }

  getSolver()->projectionBasisIs(*projectionBasis_);
  getSolver()->factor(); // Delayed factorization
}

void
GaussNewtonNonLinDynamic::postProcess() {
  std::printf("GaussNewtonNonLinDynamic::postProcess()\n");
}

GalerkinProjectionSolver *
GaussNewtonNonLinDynamic::getSolver() {
  return static_cast<GalerkinProjectionSolver *>(NonLinDynamic::getSolver());
}

bool
GaussNewtonNonLinDynamic::factorWhenBuilding() const {
  return false; // Delayed factorization
}
