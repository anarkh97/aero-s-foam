#include "DynamStateReductorImpl.h"

#include <Math.d/SymFullMatrix.h>

namespace Pita {

DynamStateReductorImpl::DynamStateReductorImpl(const DynamStateBasis * reductionBasis,
                                               const RankDeficientSolver * solver) :
  reductionBasis_(reductionBasis),
  solver_(solver)
{
  reset();
}

void
DynamStateReductorImpl::reset() {
  size_t newVectorSize = this->reductionBasis()->vectorSize();
  size_t newReducedBasisSize = this->solver()->matrixSize(); 

  if ((newVectorSize != this->vectorSize()) || (newReducedBasisSize != this->reducedBasisSize())) {
    setReducedBasisSize(newReducedBasisSize);
    setInitialState(DynamState());
  }

  setVectorSize(newVectorSize);
}
 
void
DynamStateReductorImpl::initialStateIs(const DynamState & is) {
  if (is.vectorSize() == 0) {
    getReducedBasisComponents().zero();
    setInitialState(DynamState());
    return;
  } 
  
  if (is.vectorSize() != this->vectorSize()) {
    throw Fwk::RangeException("in SimpleDynamStateReductor::initialStateIs"); 
  }

  //log() << "DynamStateReductorImpl::initialStateIs() with factorRank / reducedBasisSize = " << solver()->factorRank() << " / " << reducedBasisSize() << "\n";

  if (this->reducedBasisSize() > 0) {
    // Assemble rhs
    for (int i = 0; i < reducedBasisSize(); ++i) {
      this->getReducedBasisComponents()[i] = is * this->reductionBasis()->state(i);
    }

    // Perform in place resolution
    this->solver()->solution(this->getReducedBasisComponents());
  }

  setInitialState(is);
}

DynamStateReductorImpl::Manager::Manager(const DynamStateBasis * defaultReductionBasis,
                                         const RankDeficientSolver * defaultSolver) :
  instance_(),
  defaultReductionBasis_(defaultReductionBasis),
  defaultSolver_(defaultSolver)
{}


DynamStateReductorImpl *
DynamStateReductorImpl::Manager::instance(const String & key) const {
  InstanceContainer::const_iterator it = instance_.find(key);
  return (it != instance_.end()) ? it->second.ptr() : NULL;
}

size_t
DynamStateReductorImpl::Manager::instanceCount() const {
  return instance_.size();
}

DynamStateReductorImpl *
DynamStateReductorImpl::Manager::instanceNew(const String & key) {
  InstanceContainer::iterator it = instance_.lower_bound(key);
  if (it != instance_.end() && it->first == key)
    throw NameInUseException("in DynamStateReductorImpl::instanceNew");

  DynamStateReductorImpl::Ptr newInstance =
      new DynamStateReductorImpl(this->defaultReductionBasis(), this->defaultSolver());
  instance_.insert(it, std::make_pair(key, newInstance));
  
  return newInstance.ptr();
}

void
DynamStateReductorImpl::Manager::instanceDel(const String & key) {
  instance_.erase(key);
}

/*void
DynamStateReductorImpl::Manager::normalMatrixIs(const SymFullMatrix & nm) {
  int newSize = nm.dim();
  
  if (newSize > 0) {
    // Update solver
    solver_->matrixIs(nm);
    solver_->statusIs(RankDeficientSolver::FACTORIZED);
  }

  // Update Reductor instances 
  for (InstanceContainer::iterator it = instance_.begin(); it != instance_.end(); ++it) {
    DynamStateReductorImpl * reductor = it->second.ptr();
    
    reductor->setReducedBasisSize(newSize);
    reductor->getReducedBasisComponents().initialize(newSize);
    reductor->resetInitialState();
  }

  normalMatrix_ = &nm;
}*/

void
DynamStateReductorImpl::Manager::defaultReductionBasisIs(const DynamStateBasis * mb) {
  defaultReductionBasis_ = mb;
  resetInstances();
}

void
DynamStateReductorImpl::Manager::defaultSolverIs(const RankDeficientSolver * s) {
  defaultSolver_ = s;
  resetInstances();
}

void
DynamStateReductorImpl::Manager::resetInstances() {
  for (InstanceContainer::iterator it = instance_.begin(); it != instance_.end(); ++it) {
    it->second->reset();
  }
}

/*const SymFullMatrix *
DynamStateReductorImpl::Manager::defaultNormalMatrix() {
  static SymFullMatrix defaultNormalMatrix_;
  return &defaultNormalMatrix_;
}*/

} // end namespace Pita
