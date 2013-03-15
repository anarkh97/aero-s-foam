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

  virtual void lastMidTimeIs(double) = 0;
  virtual void lastDeltaIs(double) = 0;
  virtual void stateSnapshotAdd(const GeomState &) = 0;
  virtual void velocSnapshotAdd(const Vector &) = 0;
  virtual void accelSnapshotAdd(const Vector &) = 0;
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
  class sttSnapImpl;
  class velSnapImpl;
  class accSnapImpl;
  class resSnapImpl;
  class jacSnapImpl;

private:
  // Dummy constructor
  PodProjectionNonLinDynamicDetail();
};

// Basic implementation
class PodProjectionNonLinDynamicDetail::BasicImpl : public PodProjectionNonLinDynamic::Impl {
public:
  explicit BasicImpl(PodProjectionNonLinDynamic *parent);
  
  // Overriden functions
  virtual void lastMidTimeIs(double t);
  virtual void lastDeltaIs(double dt);
  virtual void stateSnapshotAdd(const GeomState &);
  virtual void velocSnapshotAdd(const Vector &);
  virtual void accelSnapshotAdd(const Vector &);
  virtual void handleResidualSnapshot(const Vector &res);
  virtual void handleJacobianSnapshot();

  ~BasicImpl();


protected: 
  VecNodeDof6Conversion vecNodeDof6Conversion_;
  FileNameInfo fileInfo_;

  VecBasis projectionBasis_;
};

PodProjectionNonLinDynamicDetail::BasicImpl::~BasicImpl() {}

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
  
  filePrint(stderr, " ... Projection subspace of dimension = %d ...\n", projectionBasis_.vectorCount());

  // Setup solver
  PodProjectionSolver *solver = getSolver();
  solver->projectionBasisIs(projectionBasis_);
  solver->factor(); // Delayed factorization
}

void
PodProjectionNonLinDynamicDetail::BasicImpl::lastMidTimeIs(double t) {
  //empty
}

void
PodProjectionNonLinDynamicDetail::BasicImpl::lastDeltaIs(double dt) {
  //empty
}

void
PodProjectionNonLinDynamicDetail::BasicImpl::stateSnapshotAdd(const GeomState &snap) {
 //empty
}

void
PodProjectionNonLinDynamicDetail::BasicImpl::velocSnapshotAdd(const Vector &) {
  // Nothing to do
}

void
PodProjectionNonLinDynamicDetail::BasicImpl::accelSnapshotAdd(const Vector &) {
  // Nothing to do
}

void
PodProjectionNonLinDynamicDetail::BasicImpl::handleResidualSnapshot(const Vector &) {
  // Nothing to do
}

void
PodProjectionNonLinDynamicDetail::BasicImpl::handleJacobianSnapshot() {
  // Nothing to do
}

class PodProjectionNonLinDynamicDetail::sttSnapImpl : public PodProjectionNonLinDynamicDetail::BasicImpl {
  public:
     void lastMidTimeIs(double t);
     void lastDeltaIs(double dt);
     void stateSnapshotAdd(const GeomState &state);
     void velocSnapshotAdd(const Vector &res);
     void accelSnapshotAdd(const Vector &res);
     void handleResidualSnapshot(const Vector &res);
     void handleJacobianSnapshot(); 

    int dofSetNodeCount() const { return converter_.dofSetNodeCount(); }
    int vectorSize() const { return converter_.vectorSize(); }

    explicit sttSnapImpl(Domain * domain, PodProjectionNonLinDynamic * parent);


  protected:

    Domain * domain_;
   
    template <typename VecType>
    void fillSnapBuffer(const VecType &origin);

    const NodeDof6Buffer &snapBuffer() const { return snapBuffer_; }
    const VecNodeDof6Conversion &converter() const { return converter_; }

  private:

    VecNodeDof6Conversion converter_;
    NodeDof6Buffer snapBuffer_;

    double timeStamp_;

  protected:
    BasisBinaryOutputFile stateSnapFile_;
};

// Implementation with velocity snapshots
class PodProjectionNonLinDynamicDetail::velSnapImpl : public PodProjectionNonLinDynamicDetail::BasicImpl {
public:
  explicit velSnapImpl(PodProjectionNonLinDynamic *parent);

  // Overriden functions
   void lastMidTimeIs(double t);
   void lastDeltaIs(double dt);
   void stateSnapshotAdd(const GeomState &state);
   void velocSnapshotAdd(const Vector &res);
   void accelSnapshotAdd(const Vector &res);
   void handleResidualSnapshot(const Vector &res);
   void handleJacobianSnapshot();

private:
  BasisOutputStream velocitySnapFile_;
};

// Implementation with acceleration snapshots
class PodProjectionNonLinDynamicDetail::accSnapImpl : public PodProjectionNonLinDynamicDetail::BasicImpl {
public:
  explicit accSnapImpl(PodProjectionNonLinDynamic *parent);

  // Overriden functions
   void lastMidTimeIs(double t);
   void lastDeltaIs(double dt);
   void stateSnapshotAdd(const GeomState &state);
   void velocSnapshotAdd(const Vector &res);
   void accelSnapshotAdd(const Vector &res);
   void handleResidualSnapshot(const Vector &res);
   void handleJacobianSnapshot();

private:
  BasisOutputStream accelerationSnapFile_;
};

// Implementation with residual snapshots
class PodProjectionNonLinDynamicDetail::resSnapImpl : public PodProjectionNonLinDynamicDetail::BasicImpl {
public:
  explicit resSnapImpl(PodProjectionNonLinDynamic *parent);
  
  // Overriden functions
   void lastMidTimeIs(double t);
   void lastDeltaIs(double dt);
   void stateSnapshotAdd(const GeomState &state);
   void velocSnapshotAdd(const Vector &res);
   void accelSnapshotAdd(const Vector &res);
   void handleResidualSnapshot(const Vector &res);
   void handleJacobianSnapshot();

private:
  BasisOutputStream residualSnapFile_;
};

//Implementation with jacobian snapshots
class PodProjectionNonLinDynamicDetail::jacSnapImpl : public PodProjectionNonLinDynamicDetail::BasicImpl {
public:
  explicit jacSnapImpl(PodProjectionNonLinDynamic *parent);

  // Overriden functions
   void lastMidTimeIs(double t);
   void lastDeltaIs(double dt);
   void stateSnapshotAdd(const GeomState &state);
   void velocSnapshotAdd(const Vector &res);
   void accelSnapshotAdd(const Vector &res);
   void handleResidualSnapshot(const Vector &res);
   void handleJacobianSnapshot();

private:
  BasisOutputStream jacobianSnapFile_;
};

PodProjectionNonLinDynamicDetail::sttSnapImpl::sttSnapImpl(Domain * domain, PodProjectionNonLinDynamic * parent) :
  PodProjectionNonLinDynamicDetail::BasicImpl(parent),
  domain_(domain),
  converter_(*domain->getCDSA()),
  snapBuffer_(dofSetNodeCount()),
  stateSnapFile_(BasisFileId(fileInfo_, BasisId::STATE, BasisId::SNAPSHOTS), dofSetNodeCount()),
  timeStamp_(domain->solInfo().initialTime)
{}

template <typename VecType>
inline
void
PodProjectionNonLinDynamicDetail::sttSnapImpl::fillSnapBuffer(const VecType &snap) {
  converter_.paddedNodeDof6(snap, snapBuffer_);
}

PodProjectionNonLinDynamicDetail::velSnapImpl::velSnapImpl(PodProjectionNonLinDynamic *parent) :
  PodProjectionNonLinDynamicDetail::BasicImpl(parent),
  velocitySnapFile_(BasisFileId(fileInfo_, BasisId::VELOCITY, BasisId::SNAPSHOTS), vecNodeDof6Conversion_)
{}

PodProjectionNonLinDynamicDetail::accSnapImpl::accSnapImpl(PodProjectionNonLinDynamic *parent) :
  PodProjectionNonLinDynamicDetail::BasicImpl(parent),
  accelerationSnapFile_(BasisFileId(fileInfo_, BasisId::ACCELERATION, BasisId::SNAPSHOTS), vecNodeDof6Conversion_)
{}

PodProjectionNonLinDynamicDetail::resSnapImpl::resSnapImpl(PodProjectionNonLinDynamic *parent) :
  PodProjectionNonLinDynamicDetail::BasicImpl(parent),
  residualSnapFile_(BasisFileId(fileInfo_, BasisId::RESIDUAL, BasisId::SNAPSHOTS), vecNodeDof6Conversion_)
{}

PodProjectionNonLinDynamicDetail::jacSnapImpl::jacSnapImpl(PodProjectionNonLinDynamic *parent) :
  PodProjectionNonLinDynamicDetail::BasicImpl(parent),
  jacobianSnapFile_(BasisFileId(fileInfo_, BasisId::JACOBIAN, BasisId::SNAPSHOTS), vecNodeDof6Conversion_)
{}

void
PodProjectionNonLinDynamicDetail::sttSnapImpl::lastMidTimeIs(double t) {
  timeStamp_ = t;
}

void
PodProjectionNonLinDynamicDetail::sttSnapImpl::lastDeltaIs(double dt) {
  timeStamp_ += dt;
}

void
PodProjectionNonLinDynamicDetail::velSnapImpl::lastMidTimeIs(double t) {
  //empty
}

void
PodProjectionNonLinDynamicDetail::velSnapImpl::lastDeltaIs(double dt) {
  //empty
}

void
PodProjectionNonLinDynamicDetail::accSnapImpl::lastMidTimeIs(double t) {
  //empty
}

void
PodProjectionNonLinDynamicDetail::accSnapImpl::lastDeltaIs(double dt) {
  //empty
}

void
PodProjectionNonLinDynamicDetail::resSnapImpl::lastMidTimeIs(double t) {
  //empty
}

void
PodProjectionNonLinDynamicDetail::resSnapImpl::lastDeltaIs(double dt) {
  //empty
}

void
PodProjectionNonLinDynamicDetail::jacSnapImpl::lastMidTimeIs(double t) {
  //empty
}

void
PodProjectionNonLinDynamicDetail::jacSnapImpl::lastDeltaIs(double dt) {
  //empty
}

void
PodProjectionNonLinDynamicDetail::sttSnapImpl::stateSnapshotAdd(const GeomState &snap) {
  const CoordSet &refCoords = domain_->getNodes();

  for (int iNode = 0, iNodeEnd = dofSetNodeCount(); iNode != iNodeEnd; ++iNode) {
    double *nodeBuffer = snapBuffer_[iNode];

    const Node *refNode = refCoords[iNode];
    if (refNode) {
      const NodeState &snapNode = snap[iNode];

      // Translational dofs
      nodeBuffer[0] = snapNode.x - refNode->x;
      nodeBuffer[1] = snapNode.y - refNode->y;
      nodeBuffer[2] = snapNode.z - refNode->z;

      // Rotational dofs
      mat_to_vec(const_cast<double (*)[3]>(snapNode.R), &nodeBuffer[3]);
    } else {
      // Node does not really exist, corresponds to a gap in node numbering
      std::fill_n(nodeBuffer, 6, 0.0);
    }
  }

  stateSnapFile_.stateAdd(snapBuffer_, timeStamp_);
}

void
PodProjectionNonLinDynamicDetail::velSnapImpl::stateSnapshotAdd(const GeomState &snap) {
 //empty
}

void
PodProjectionNonLinDynamicDetail::accSnapImpl::stateSnapshotAdd(const GeomState &snap) {
  //empty
}

void
PodProjectionNonLinDynamicDetail::resSnapImpl::stateSnapshotAdd(const GeomState &snap) {
 //empty
}

void
PodProjectionNonLinDynamicDetail::jacSnapImpl::stateSnapshotAdd(const GeomState &snap) {
  //empty
}

void
PodProjectionNonLinDynamicDetail::sttSnapImpl::velocSnapshotAdd(const Vector &veloc) {
  //empty
}

void
PodProjectionNonLinDynamicDetail::velSnapImpl::velocSnapshotAdd(const Vector &veloc) {
  velocitySnapFile_ << veloc;
}

void
PodProjectionNonLinDynamicDetail::accSnapImpl::velocSnapshotAdd(const Vector &veloc) {
  //empty
}

void
PodProjectionNonLinDynamicDetail::resSnapImpl::velocSnapshotAdd(const Vector &veloc) {
  //empty
}

void
PodProjectionNonLinDynamicDetail::jacSnapImpl::velocSnapshotAdd(const Vector &veloc) {
  //empty
}

void
PodProjectionNonLinDynamicDetail::sttSnapImpl::accelSnapshotAdd(const Vector &accel) {
  //empty
}

void
PodProjectionNonLinDynamicDetail::velSnapImpl::accelSnapshotAdd(const Vector &accel) {
 //empty
}

void
PodProjectionNonLinDynamicDetail::accSnapImpl::accelSnapshotAdd(const Vector &accel) {
  accelerationSnapFile_ << accel;
}

void
PodProjectionNonLinDynamicDetail::resSnapImpl::accelSnapshotAdd(const Vector &accel) {
  //empty
}

void
PodProjectionNonLinDynamicDetail::jacSnapImpl::accelSnapshotAdd(const Vector &accel) {
  //empty
}

void
PodProjectionNonLinDynamicDetail::sttSnapImpl::handleResidualSnapshot(const Vector &res) {
  //empty
}

void
PodProjectionNonLinDynamicDetail::velSnapImpl::handleResidualSnapshot(const Vector &res) {
  //empty
}

void
PodProjectionNonLinDynamicDetail::accSnapImpl::handleResidualSnapshot(const Vector &res) {
  //empty
}

void
PodProjectionNonLinDynamicDetail::resSnapImpl::handleResidualSnapshot(const Vector &res) {
  residualSnapFile_ << res;
}

void
PodProjectionNonLinDynamicDetail::jacSnapImpl::handleResidualSnapshot(const Vector &res) {
  //empty
}

void
PodProjectionNonLinDynamicDetail::sttSnapImpl::handleJacobianSnapshot() {
  //empty
}

void
PodProjectionNonLinDynamicDetail::velSnapImpl::handleJacobianSnapshot() {
  //empty
}

void
PodProjectionNonLinDynamicDetail::accSnapImpl::handleJacobianSnapshot() {
  //empty
}

void
PodProjectionNonLinDynamicDetail::resSnapImpl::handleJacobianSnapshot() {
  //empty
}

void
PodProjectionNonLinDynamicDetail::jacSnapImpl::handleJacobianSnapshot() {
  Vector snap(solVecInfo());
  expand(getSolver()->lastReducedMatrixAction(), getSolver()->lastReducedSolution(), snap);
  jacobianSnapFile_ << snap;
}


// Main class members

PodProjectionNonLinDynamic::PodProjectionNonLinDynamic(Domain *d) :
  NonLinDynamic(d),
  impl_(NULL),
  sttImpl_(NULL),
  velImpl_(NULL),
  accImpl_(NULL),
  resImpl_(NULL),
  jacImpl_(NULL)
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
   if(domain->solInfo().statevectPodRom)
    sttImpl_.reset(new PodProjectionNonLinDynamicDetail::sttSnapImpl(this->domain,this)); 
   if(domain->solInfo().velocvectPodRom)
    velImpl_.reset(new PodProjectionNonLinDynamicDetail::velSnapImpl(this));
   if(domain->solInfo().accelvectPodRom)
    accImpl_.reset(new PodProjectionNonLinDynamicDetail::accSnapImpl(this));
   if(domain->solInfo().residvectPodRom) 
    resImpl_.reset(new PodProjectionNonLinDynamicDetail::resSnapImpl(this));
   if(domain->solInfo().jacobvectPodRom)
    jacImpl_.reset(new PodProjectionNonLinDynamicDetail::jacSnapImpl(this));
  } else {
    impl_.reset(new PodProjectionNonLinDynamicDetail::BasicImpl(this));
  }
}


const PodProjectionSolver *
PodProjectionNonLinDynamic::getSolver() const {
  return dynamic_cast<PodProjectionSolver *>(const_cast<PodProjectionNonLinDynamic *>(this)->NonLinDynamic::getSolver());
}

PodProjectionSolver *
PodProjectionNonLinDynamic::getSolver() {
  return const_cast<PodProjectionSolver *>(const_cast<const PodProjectionNonLinDynamic *>(this)->getSolver());
}

int
PodProjectionNonLinDynamic::checkConvergence(int iteration, double normRes, Vector &residual, Vector &dv, double time) {
  if(domain->solInfo().jacobvectPodRom)
  jacImpl_->handleJacobianSnapshot();

  // Forward to hidden base class function
  return NonLinDynamic::checkConvergence(iteration, normRes, residual, dv, time); 
}

double
PodProjectionNonLinDynamic::getResidualNorm(const Vector &residual, GeomState &, double) {
  return getSolver()->projectAndComputeNorm(residual);
}

void
PodProjectionNonLinDynamic::handleResidualSnapshot(const Vector &snap) {
  if(domain->solInfo().residvectPodRom)
  resImpl_->handleResidualSnapshot(snap);
}

bool
PodProjectionNonLinDynamic::factorWhenBuilding() const {
  return false; // Delayed factorization
}

void
PodProjectionNonLinDynamic::saveMidTime(double t) {
  if(domain->solInfo().statevectPodRom)
  sttImpl_->lastMidTimeIs(t);
}

void
PodProjectionNonLinDynamic::saveDelta(double dt) {
  if(domain->solInfo().statevectPodRom)
  sttImpl_->lastDeltaIs(dt);
}

void
PodProjectionNonLinDynamic::saveStateSnapshot(const GeomState &state) {
  if(domain->solInfo().statevectPodRom)
  sttImpl_->stateSnapshotAdd(state);
}

void
PodProjectionNonLinDynamic::saveVelocitySnapshot(const Vector &veloc) {
  if(domain->solInfo().velocvectPodRom)
  velImpl_->velocSnapshotAdd(veloc);
}

void
PodProjectionNonLinDynamic::saveAccelerationSnapshot(const Vector &accel) {
  if(domain->solInfo().accelvectPodRom)
  accImpl_->accelSnapshotAdd(accel);
}

} /* end namespace Rom */
