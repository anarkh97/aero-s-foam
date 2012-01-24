#include "PodProjectionNonLinDynamic.h"

#include "FileNameInfo.h"
#include "BasisFileStream.h"
#include "NodeDof6Buffer.h"
#include "VecNodeDof6Conversion.h"

#include "VecBasis.h"
#include "VecBasisFile.h"

#include "BasisOps.h"

#include <Driver.d/Domain.h>
#include <Utils.d/DistHelper.h>

#include <algorithm>
#include <stdexcept>
#include <cstddef>

namespace Rom {

// Implementation classes

// Common implementation
class PodProjectionNonLinDynamic::Impl {
public:
  virtual void handleResidualSnapshot(const Vector &res) = 0;
  virtual void handleJacobianSnapshot() = 0;

  virtual ~Impl();

protected:
  explicit Impl(PodProjectionNonLinDynamic *parent);

  const ConstrainedDSA &getCDSA() const { return *parent_->domain->getCDSA(); }
  int solVecInfo() const { return parent_->solVecInfo(); }
  const SolverInfo &solInfo() const { return parent_->domain->solInfo(); }
  PodProjectionSolver *getSolver() { return parent_->getSolver(); }
  const PodProjectionSolver *getSolver() const { return parent_->getSolver(); }
    
private:  
  PodProjectionNonLinDynamic *parent_;
  
  // Disallow copy & assignment
  Impl(const Impl &);
  Impl &operator=(const Impl &);
};

PodProjectionNonLinDynamic::Impl::Impl(PodProjectionNonLinDynamic *parent) :
  parent_(parent)
{}

PodProjectionNonLinDynamic::Impl::~Impl() {
  // Nothing to do
}

// Dummy class, used for namespace access
class PodProjectionNonLinDynamicDetail : private PodProjectionNonLinDynamic {
public:
  class BasicImpl;
  class SnapImpl;

private:
  // Dummy constructor
  PodProjectionNonLinDynamicDetail();
};

// Basic implementation
class PodProjectionNonLinDynamicDetail::BasicImpl : public PodProjectionNonLinDynamic::Impl {
public:
  explicit BasicImpl(PodProjectionNonLinDynamic *parent);
  
  // Overriden functions
  virtual void handleResidualSnapshot(const Vector &res);
  virtual void handleJacobianSnapshot();

protected: 
  VecNodeDof6Conversion vecNodeDof6Conversion_;
  FileNameInfo fileInfo_;

  VecBasis projectionBasis_;
};

PodProjectionNonLinDynamicDetail::BasicImpl::BasicImpl(PodProjectionNonLinDynamic *parent) :
  PodProjectionNonLinDynamic::Impl(parent),
  vecNodeDof6Conversion_(getCDSA()),
  fileInfo_()
{
  // Load projection basis
  BasisInputStream projectionBasisInput(BasisFileId(fileInfo_, BasisId::STATE, BasisId::POD), vecNodeDof6Conversion_);

  if (projectionBasisInput.vectorSize() != solVecInfo()) {
    throw std::domain_error("Projection basis has incorrect #rows");
  }

  const int projectionSubspaceSize = solInfo().maxSizePodRom ?
                                     std::min(solInfo().maxSizePodRom, projectionBasisInput.size()) :
                                     projectionBasisInput.size();
  
  readVectors(projectionBasisInput, projectionBasis_, projectionSubspaceSize);
  
  filePrint(stderr, "Projection subspace of dimension = %d\n", projectionBasis_.vectorCount());

  // Setup solver
  PodProjectionSolver *solver = getSolver();
  solver->projectionBasisIs(projectionBasis_);
  solver->factor(); // Delayed factorization
}

void
PodProjectionNonLinDynamicDetail::BasicImpl::handleResidualSnapshot(const Vector &) {
  // Nothing to do
}

void
PodProjectionNonLinDynamicDetail::BasicImpl::handleJacobianSnapshot() {
  // Nothing to do
}

// Implementation with residual/jacobian snapshots
class PodProjectionNonLinDynamicDetail::SnapImpl : public PodProjectionNonLinDynamicDetail::BasicImpl {
public:
  explicit SnapImpl(PodProjectionNonLinDynamic *parent);
  
  // Overriden functions
  virtual void handleResidualSnapshot(const Vector &res);
  virtual void handleJacobianSnapshot();

private:
  BasisOutputStream residualSnapFile_;
  BasisOutputStream jacobianSnapFile_;
};

PodProjectionNonLinDynamicDetail::SnapImpl::SnapImpl(PodProjectionNonLinDynamic *parent) :
  PodProjectionNonLinDynamicDetail::BasicImpl(parent),
  residualSnapFile_(BasisFileId(fileInfo_, BasisId::RESIDUAL, BasisId::SNAPSHOTS), vecNodeDof6Conversion_),
  jacobianSnapFile_(BasisFileId(fileInfo_, BasisId::JACOBIAN, BasisId::SNAPSHOTS), vecNodeDof6Conversion_)
{}

void
PodProjectionNonLinDynamicDetail::SnapImpl::handleResidualSnapshot(const Vector &res) {
  residualSnapFile_ << res;
}

void
PodProjectionNonLinDynamicDetail::SnapImpl::handleJacobianSnapshot() {
  Vector snap(solVecInfo());
  expand(getSolver()->lastReducedMatrixAction(), getSolver()->lastReducedSolution(), snap);
  jacobianSnapFile_ << snap;
}


// Main class members

PodProjectionNonLinDynamic::PodProjectionNonLinDynamic(Domain *d) :
  NonLinDynamic(d),
  impl_(NULL)
{}

PodProjectionNonLinDynamic::~PodProjectionNonLinDynamic() {
  // Nothing to do
}

void
PodProjectionNonLinDynamic::preProcess() {
  NonLinDynamic::preProcess();
  
  if (!dynamic_cast<PodProjectionSolver *>(NonLinDynamic::getSolver())) {
    throw std::runtime_error("Solver must be a PodProjectionSolver");
  }

  if (domain->solInfo().snapshotsPodRom) {
    impl_.reset(new PodProjectionNonLinDynamicDetail::SnapImpl(this));
  } else {
    impl_.reset(new PodProjectionNonLinDynamicDetail::BasicImpl(this));
  }
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
  impl_->handleJacobianSnapshot();

  // Forward to hidden base class function
  return NonLinDynamic::checkConvergence(iteration, normRes, residual, dv, time); 
}

double
PodProjectionNonLinDynamic::getResidualNorm(const Vector &residual) {
  return getSolver()->projectAndComputeNorm(residual);
}

void
PodProjectionNonLinDynamic::handleResidualSnapshot(const Vector &snap) {
  impl_->handleResidualSnapshot(snap);
}

bool
PodProjectionNonLinDynamic::factorWhenBuilding() const {
  return false; // Delayed factorization
}

} /* end namespace Rom */
