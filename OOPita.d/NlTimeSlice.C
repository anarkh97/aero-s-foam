#include "NlTimeSlice.h"
#include "PitaNonLinDynam.h"
#include "NlDynamTimeIntegrator.h"
#include "AlternatingIntegratorPropagator.h"
#include "LinearizedPropagator.h"
#include "NlDynamOps.h"

namespace Pita {

//
// NlTimeSlice implementation
//

// Constructor

NlTimeSlice::NlTimeSlice(SliceRank rank, Seconds initialTime, Seconds finalTime, Manager * manager, Status status) :
  LocalTimeSlice(rank, initialTime, finalTime, status), manager_(manager) {}

// Public members

void
NlTimeSlice::initialProjectionBasisInc(DynamStateBasis::PtrConst basis) {
  if (status() == TimeSlice::inactive || status() == active) {
    seedUpdatePropagator()->projectionBasisInc(basis);
    setInitialProjectionBasis(seedUpdatePropagator()->projectionBasis().ptr());
  }
}

// Protected members

void
NlTimeSlice::activateSlice() {
  seedUpdatePropagator()->referenceDisplacementIs(initialInterface()->seed().displacement(), initialTime());
}

void
NlTimeSlice::optimizeLastIteration() {
  seedUpdatePropagator()->projectionBasisIs(DynamStatePlainBasis::New(seedUpdatePropagator()->projectionBasis()->vectorSize()));
}

//
// Manager implementation
//

// Public members  

NlTimeSlice::Manager::Manager(NlDynamTimeIntegrator * integrator, CommManager * commManager) :
  LocalTimeSlice::Manager(commManager),
  probDesc_(const_cast<PitaNonLinDynamic*>(integrator->probDesc())),
  dynamOps_(NlDynamOps::New(probDesc_).ptr()),
  timeIntegrator_(integrator)
{
  setOutgoingInitialStateBasis(DynamStatePlainBasis::New(probDesc_->solVecInfo()).ptr());
}

void
NlTimeSlice::Manager::timeSliceDel(SliceRank rank) {
  SliceMap::iterator it = slice_.find(rank);
  if (it != slice_.end())
    slice_.erase(it);
}

// Protected members

NlTimeSlice *
NlTimeSlice::Manager::createTimeSlice(SliceRank rank) {
  // Find insertion point
  SliceMap::iterator it = slice_.lower_bound(rank);
  if (it != slice_.end() && it->first == rank)
    throw Fwk::NameInUseException();
  
  // Build timeslice
  Seconds t0 = Seconds(rank.value() * coarseTimeStep().value());
  Seconds tf = t0 + coarseTimeStep(); 
  NlTimeSlice::Ptr newSlice = NlTimeSlice::New(rank, t0, tf, this);
  
  // Set propagators
  AlternatingIntegratorPropagator::Ptr finePropagator = AlternatingIntegratorPropagator::New(timeIntegrator_.ptr());
  finePropagator->initialTimeIs(t0);
  finePropagator->timeStepCountIs(timeGridRatio());
  finePropagator->sliceRankIs(rank);
  newSlice->setLocalPropagator(finePropagator.ptr());
  DynamPropagator::Ptr linearizedFinePropagator = LinearizedPropagator::New(timeIntegrator_.ptr());
  ProjectorPropagator::Ptr projectorPropagator = ProjectorPropagator::New(linearizedFinePropagator.ptr(), dynamOps_.ptr());
  finePropagator->secondaryPropagatorIs(projectorPropagator);
  newSlice->setSeedUpdatePropagator(projectorPropagator.ptr());

  // Set projection base
  newSlice->setInitialProjectionBasis(DynamStatePlainBasis::New(vectorSize()).ptr());
  
  // Add timeslice 
  slice_.insert(it, std::make_pair(rank, newSlice));
  return newSlice.ptr(); 
}

NlTimeSlice *
NlTimeSlice::Manager::findTimeSlice(SliceRank rank) const {
  SliceMap::const_iterator it = slice_.find(rank);
  return (it != slice_.end()) ? it->second.ptr() : NULL;
}

NlTimeSlice *
NlTimeSlice::Manager::getNextTimeSlice(const TimeSlice * ts) const {
  SliceMap::const_iterator it = ts ? slice_.upper_bound(ts->rank()) : slice_.begin();
  return (it != slice_.end()) ? it->second.ptr() : NULL;
}

TimeStepCount
NlTimeSlice::Manager::timeGridRatio() const {
  return TimeStepCount(probDesc()->getJratio());
}

Seconds
NlTimeSlice::Manager::coarseTimeStep() const {
  return Seconds(probDesc()->getCoarseDt());
}

size_t
NlTimeSlice::Manager::vectorSize() const {
  return probDesc()->solVecInfo();
}

} // end namespace Pita
