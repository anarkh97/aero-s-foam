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
  
  if (!dynamic_cast<GaussNewtonSolver *>(NonLinDynamic::getSolver())) {
    throw std::runtime_error("Solver must be a GaussNewtonSolver");
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


const GaussNewtonSolver *
GaussNewtonNonLinDynamic::getSolver() const {
  return static_cast<GaussNewtonSolver *>(const_cast<GaussNewtonNonLinDynamic *>(this)->NonLinDynamic::getSolver());
}

GaussNewtonSolver *
GaussNewtonNonLinDynamic::getSolver() {
  return const_cast<GaussNewtonSolver *>(const_cast<const GaussNewtonNonLinDynamic *>(this)->getSolver());
}

int
GaussNewtonNonLinDynamic::checkConvergence(int iteration, double normRes, Vector &residual, Vector &dv, double time) {
  computeAndSaveJacobianSnapshot();

  // Forward to hidden base class function
  return NonLinDynamic::checkConvergence(iteration, normRes, residual, dv, time); 
}

double
GaussNewtonNonLinDynamic::getResidualNorm(const Vector &residual) {
  return getSolver()->projectAndComputeNorm(residual);
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
