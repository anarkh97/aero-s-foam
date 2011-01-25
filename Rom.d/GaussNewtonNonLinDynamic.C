#include "GaussNewtonNonLinDynamic.h"

#include "BasisFileIterator.h"

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
  BasisInputRange projectionBasisInput("GaussNewtonBasis", *vecNodeDof6Conversion_); //TODO: file name
  assert(projectionBasisInput.vectorSize() == solVecInfo());
  std::fprintf(stderr, "Gauss-Newton projection basis size = %d\n", projectionBasisInput.size());

  projectionBasis_.reset(new GenVecBasis<double>(projectionBasisInput.size(), projectionBasisInput.vectorSize()));

  GenVecBasis<double>::iterator jt = projectionBasis_->begin();
  for (BasisInputRange::const_iterator it = projectionBasisInput.begin();
      it != projectionBasisInput.end(); ++it) {
    (*it)(*jt++);
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
  Vector jacobianSnap(solVecInfo());
  expand(getSolver()->lastReducedMatrixAction(), getSolver()->lastReducedSolution(), jacobianSnap);
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
