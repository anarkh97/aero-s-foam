#include "GaussNewtonNonLinDynamic.h"

#include <Driver.d/Domain.h>

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

  vecNodeDof6Conversion_.reset(new VecNodeDof6Conversion(*this->domain->getCDSA()));
  snapBuffer_.sizeIs(vecNodeDof6Conversion_->nodeCount());

  residualSnapFile_.reset(new BasisOutputFile("RomResidualSnap", vecNodeDof6Conversion_->nodeCount())); //TODO name
  jacobianSnapFile_.reset(new BasisOutputFile("RomJacobianSnap", vecNodeDof6Conversion_->nodeCount())); //TODO name
}

void
GaussNewtonNonLinDynamic::postProcess() {
  residualSnapFile_->updateStateCountStatus();
  jacobianSnapFile_->updateStateCountStatus();
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
  const GalerkinProjectionSolver *solver = getSolver();

  Vector jacobianSnap(solVecInfo(), 0.0);
  const int iEnd = solver->basisSize();
  for (int i = 0; i < iEnd; ++i) {
    jacobianSnap.linAdd(solver->lastReducedSolution()[i], solver->lastReducedMatrixAction()[i]);
  }

  saveJacobianSnapshot(jacobianSnap);
}

void
GaussNewtonNonLinDynamic::saveResidualSnapshot(const Vector &snap) {
  residualSnapFile_->stateAdd(vecNodeDof6Conversion_->nodeDof6(snap, snapBuffer_));
}

void
GaussNewtonNonLinDynamic::saveJacobianSnapshot(const Vector &snap) {
  jacobianSnapFile_->stateAdd(vecNodeDof6Conversion_->nodeDof6(snap, snapBuffer_));
}

bool
GaussNewtonNonLinDynamic::factorWhenBuilding() const {
  return false; // Delayed factorization
}
