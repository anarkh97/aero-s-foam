#include "GaussNewtonNonLinDynamic.h"

#include "BasisFileStream.h"
#include "VecBasisFile.h"
#include "FileNameInfo.h"
#include "BasisOps.h"

#include <Driver.d/Domain.h>
#include <Utils.d/DistHelper.h>

#include <algorithm>
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

  // Input/output state conversion
  vecNodeDof6Conversion_.reset(new VecNodeDof6Conversion(*this->domain->getCDSA()));
  snapBuffer_.sizeIs(vecNodeDof6Conversion_->nodeCount());

  FileNameInfo fileInfo; 

  // Load projection basis
  BasisInputStream projectionBasisInput(fileInfo.fileName(BasisId(BasisId::STATE, BasisId::POD)), *vecNodeDof6Conversion_);

  if (projectionBasisInput.vectorSize() != solVecInfo()) {
    throw std::domain_error("Projection basis has incorrect #rows");
  }

  const int projectionSubspaceSize = domain->solInfo().maxSizePodRom ?
                                     std::min(domain->solInfo().maxSizePodRom, projectionBasisInput.size()) :
                                     projectionBasisInput.size();
  
  readVectors(projectionBasisInput, projectionBasis_, projectionSubspaceSize);
  
  filePrint(stderr, "Projection subspace of dimension = %d\n", projectionBasis_.vectorCount());

  // Setup solver
  getSolver()->projectionBasisIs(projectionBasis_);
  getSolver()->factor(); // Delayed factorization
  
  // Snapshot output
  residualSnapFile_.reset(new BasisOutputStream(fileInfo.fileName(BasisId(BasisId::RESIDUAL, BasisId::SNAPSHOTS)), *vecNodeDof6Conversion_));
  jacobianSnapFile_.reset(new BasisOutputStream(fileInfo.fileName(BasisId(BasisId::JACOBIAN, BasisId::SNAPSHOTS)), *vecNodeDof6Conversion_));
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

double
GaussNewtonNonLinDynamic::getResidualNorm(const Vector &residual) const {
  Vector reducedResidual(projectionBasis_.numVec());
  reduce(projectionBasis_, residual, reducedResidual);
  return reducedResidual.norm();
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
