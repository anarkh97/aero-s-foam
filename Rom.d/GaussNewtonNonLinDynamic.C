#include "GaussNewtonNonLinDynamic.h"

#include "BasisFileStream.h"
#include "VecBasisFile.h"
#include "BasisOps.h"

#include <Driver.d/Domain.h>

#include <stdexcept>

#include <cstdio>

GaussNewtonNonLinDynamic::GaussNewtonNonLinDynamic(Domain *d) :
  NonLinDynamic(d)
{}

void
GaussNewtonNonLinDynamic::preProcess() {
  NonLinDynamic::preProcess();
  
  if (!dynamic_cast<GalerkinProjectionSolver *>(NonLinDynamic::getSolver())) {
    throw std::runtime_error("Solver must be a GalerkinProjectionSolver");
  }

  // Input/output state conversion
  vecNodeDof6Conversion_.reset(new VecNodeDof6Conversion(*this->domain->getCDSA()));
  snapBuffer_.sizeIs(vecNodeDof6Conversion_->nodeCount());

  // Load projection basis
  BasisInputStream projectionBasisInput("GaussNewtonBasis", *vecNodeDof6Conversion_); //TODO: file name
  assert(projectionBasisInput.vectorSize() == solVecInfo());
  std::fprintf(stderr, "Gauss-Newton projection basis size = %d\n", projectionBasisInput.size());
  projectionBasisInput >> projectionBasis_;

  // Setup solver
  getSolver()->projectionBasisIs(projectionBasis_);
  getSolver()->factor(); // Delayed factorization
  
  // Snapshot output
  residualSnapFile_.reset(new BasisOutputStream("RomResidualSnap", *vecNodeDof6Conversion_)); //TODO file name
  jacobianSnapFile_.reset(new BasisOutputStream("RomJacobianSnap", *vecNodeDof6Conversion_)); //TODO file name
}

GalerkinProjectionSolver *
GaussNewtonNonLinDynamic::getSolver() {
  return static_cast<GalerkinProjectionSolver *>(NonLinDynamic::getSolver());
}

int
GaussNewtonNonLinDynamic::checkConvergence(int iteration, double normRes, Vector &residual, Vector &dv, double time) {
  computeAndSaveJacobianSnapshot();

  // Forward to hidden base class function
  return NonLinDynamic::checkConvergence(iteration, normRes, residual, dv, time); 
}

void
GaussNewtonNonLinDynamic::computeAndSaveJacobianSnapshot() {
  Vector snap(solVecInfo());
  expand(getSolver()->lastReducedMatrixAction(), getSolver()->lastReducedSolution(), snap);
  *jacobianSnapFile_ << snap;
}

bool
GaussNewtonNonLinDynamic::factorWhenBuilding() const {
  return false; // Delayed factorization
}
