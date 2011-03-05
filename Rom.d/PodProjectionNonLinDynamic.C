#include "PodProjectionNonLinDynamic.h"

#include "BasisFileStream.h"
#include "VecBasisFile.h"
#include "FileNameInfo.h"
#include "BasisOps.h"

#include <Driver.d/Domain.h>
#include <Utils.d/DistHelper.h>

#include <algorithm>
#include <stdexcept>

namespace Rom {

PodProjectionNonLinDynamic::PodProjectionNonLinDynamic(Domain *d) :
  NonLinDynamic(d)
{}

void
PodProjectionNonLinDynamic::preProcess() {
  NonLinDynamic::preProcess();
  
  if (!dynamic_cast<PodProjectionSolver *>(NonLinDynamic::getSolver())) {
    throw std::runtime_error("Solver must be a PodProjectionSolver");
  }

  // Input/output state conversion
  vecNodeDof6Conversion_.reset(new VecNodeDof6Conversion(*this->domain->getCDSA()));
  snapBuffer_.sizeIs(vecNodeDof6Conversion_->nodeCount());

  FileNameInfo fileInfo; 

  // Load projection basis
  BasisInputStream projectionBasisInput(BasisFileId(fileInfo, BasisId::STATE, BasisId::POD), *vecNodeDof6Conversion_);

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
  residualSnapFile_.reset(new BasisOutputStream(BasisFileId(fileInfo, BasisId::RESIDUAL, BasisId::SNAPSHOTS), *vecNodeDof6Conversion_));
  jacobianSnapFile_.reset(new BasisOutputStream(BasisFileId(fileInfo, BasisId::JACOBIAN, BasisId::SNAPSHOTS), *vecNodeDof6Conversion_));
}


const PodProjectionSolver *
PodProjectionNonLinDynamic::getSolver() const {
  return static_cast<PodProjectionSolver *>(const_cast<PodProjectionNonLinDynamic *>(this)->NonLinDynamic::getSolver());
}

PodProjectionSolver *
PodProjectionNonLinDynamic::getSolver() {
  return const_cast<PodProjectionSolver *>(const_cast<const PodProjectionNonLinDynamic *>(this)->getSolver());
}

int
PodProjectionNonLinDynamic::checkConvergence(int iteration, double normRes, Vector &residual, Vector &dv, double time) {
  computeAndSaveJacobianSnapshot();

  // Forward to hidden base class function
  return NonLinDynamic::checkConvergence(iteration, normRes, residual, dv, time); 
}

double
PodProjectionNonLinDynamic::getResidualNorm(const Vector &residual) {
  return getSolver()->projectAndComputeNorm(residual);
}

void
PodProjectionNonLinDynamic::computeAndSaveJacobianSnapshot() {
  Vector snap(solVecInfo());
  expand(getSolver()->lastReducedMatrixAction(), getSolver()->lastReducedSolution(), snap);
  *jacobianSnapFile_ << snap;
}

bool
PodProjectionNonLinDynamic::factorWhenBuilding() const {
  return false; // Delayed factorization
}

} /* end namespace Rom */
