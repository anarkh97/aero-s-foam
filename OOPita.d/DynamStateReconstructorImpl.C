#include "DynamStateReconstructorImpl.h"

namespace Pita {

DynamStateReconstructorImpl::DynamStateReconstructorImpl(const DynamStateBasis * rb) :
  reconstructionBasis_(rb)
{
  reset();
}

void
DynamStateReconstructorImpl::reducedBasisComponentsIs(const Vector & c) {
  int rbs = static_cast<int>(this->reducedBasisSize());
  if (c.size() != rbs) {
    throw Fwk::RangeException("in SimpleDynamStateReconstructor::reducedBasisComponentsIs");
  }
  DynamState result(this->vectorSize(), 0.0);
  for (int i = 0; i < rbs; ++i) {
    if (c[i] != 0.0) {
      result.linAdd(c[i], reconstructionBasis_->state(i));
    }
  }
  setFinalState(result);
}

void
DynamStateReconstructorImpl::reset() {
  size_t newVectorSize = this->reconstructionBasis_->vectorSize();
  size_t newReducedBasisSize = this->reconstructionBasis_->stateCount();

  if ((newVectorSize != this->vectorSize()) || (newReducedBasisSize != this->reducedBasisSize())) {
    setReducedBasisSize(newReducedBasisSize);
    setFinalState(DynamState());
  }

  setVectorSize(reconstructionBasis_->vectorSize());
}

DynamStateReconstructorImpl::Manager::Manager(const DynamStateBasis * drb) :
  defaultReconstructionBasis_(drb)
{}

DynamStateReconstructorImpl *
DynamStateReconstructorImpl::Manager::createNewInstance(const String & key) {
  return new DynamStateReconstructorImpl(this->defaultReconstructionBasis());
}

void
DynamStateReconstructorImpl::Manager::defaultReconstructionBasisIs(const DynamStateBasis * drb) {
  setDefaultReconstructionBasis(drb);
  resetInstances(); // HACK ??
}

void
DynamStateReconstructorImpl::Manager::resetInstances() {
  for (Impl::InstanceMap::iterator it = instanceBegin(); it != instanceEnd(); ++it) {
    it->second->reset();
  }
}

} // end namespace Pita
