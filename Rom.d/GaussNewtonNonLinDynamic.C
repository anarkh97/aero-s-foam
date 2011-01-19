#include "GaussNewtonNonLinDynamic.h"

#include "BasisInputFile.h"

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
  BasisInputFile projectionBasisInput("GaussNewtonBasis"); //TODO: file name
  std::fprintf(stderr, "Gauss-Newton projection basis size = %d\n", projectionBasisInput.stateCount());

  projectionBasis_.reset(new GenVectorSet<double>(projectionBasisInput.stateCount(), solVecInfo()));
  for (int i = 0; i < projectionBasis_->numVec(); ++i) {
    assert(i == projectionBasisInput.currentStateIndex());
    assert(snapBuffer_.size() == projectionBasisInput.nodeCount());
    
    projectionBasisInput.currentStateBuffer(snapBuffer_);
    
    assert(projectionBasis_->size() == vecNodeDof6Conversion_->vectorSize());
    
    vecNodeDof6Conversion_->vector(snapBuffer_, (*projectionBasis_)[i]);
    projectionBasisInput.currentStateIndexInc();
  }

  // Setup solver
  getSolver()->projectionBasisIs(*projectionBasis_);
  getSolver()->factor(); // Delayed factorization
  
  // Snapshot output
  residualSnapFile_.reset(new BasisOutputFile("RomResidualSnap", vecNodeDof6Conversion_->nodeCount())); //TODO file name
  jacobianSnapFile_.reset(new BasisOutputFile("RomJacobianSnap", vecNodeDof6Conversion_->nodeCount())); //TODO file name
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
