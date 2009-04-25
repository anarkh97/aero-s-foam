#include "HalfSliceCorrectionNetworkImpl.h"

#include "NonHomogeneousBasisCollectorImpl.h"

#include "DynamStateBasisWrapper.h"

#include <Math.d/SparseMatrix.h>
#include <Comm.d/Communicator.h>

#include <stdexcept>

#include <Timers.d/GetTime.h>

namespace Pita {

HalfSliceCorrectionNetworkImpl::HalfSliceCorrectionNetworkImpl(size_t vSize,
                                                               Communicator * timeComm,
                                                               CpuRank myCpu,
                                                               const HalfSliceSchedule * schedule,
                                                               const SliceMapping * mapping,
                                                               const DynamOps * metric,
                                                               Strategy strategy) :
  vectorSize_(vSize),
  timeCommunicator_(timeComm),
  localCpu_(myCpu),
  correctionPhase_(schedule->correction()),
  schedulingPhase_(schedule->endIteration()),
  mapping_(mapping),
  metric_(metric),
  gBuffer_(),
  mBuffer_(),
  mpiParameters_(),
  metricBasis_(DynamStatePlainBasis::New(vectorSize_)),
  finalBasis_(DynamStatePlainBasis::New(vectorSize_)),
  normalMatrix_(),
  solver_(PivotedCholeskySolver::New()),
  collector_(strategy == NON_HOMOGENEOUS ?
      ptr_cast<HalfSliceBasisCollectorImpl>(NonHomogeneousBasisCollectorImpl::New()) :
      HalfSliceBasisCollectorImpl::New()),
  reductorMgr_(DynamStateReductorImpl::Manager::New(metricBasis_.ptr(), solver_.ptr())),
  reconstructorMgr_(DynamStateReconstructorImpl::Manager::New(finalBasis_.ptr())),
  projectionBuildingReactor_(new ProjectionBuildingReactor(activityManagerInstance()->activityNew("ProjectionBuilding").ptr(), this)),
  strategy_(strategy),
  schedulingReactor_(strategy == NON_HOMOGENEOUS ?
      new NonHomogeneousSchedulingReactor(activityManagerInstance()->activityNew("ScheduleProjectionBuilding").ptr(), this) :
      new SchedulingReactor(activityManagerInstance()->activityNew("ScheduleProjectionBuilding").ptr(), this)),
  globalExchangeNumbering_(NULL)
{
  projectionBuildingReactor_->notifier()->phaseIs(correctionPhase_);
  schedulingReactor_->notifier()->phaseIs(schedulingPhase_);
  schedulingReactor_->notifier()->iterationIs(schedulingReactor_->notifier()->currentIteration());
  schedulingReactor_->notifier()->statusIs(Activity::scheduled);
}


void
HalfSliceCorrectionNetworkImpl::buildProjection() {
  double tic, toc;
  tic = getTime();

  // Build buffer
  size_t stateSize = 2 * vectorSize_;
  size_t targetBufferSize = stateSize * globalExchangeNumbering_->stateCount();
  if (gBuffer_.size() < targetBufferSize) {
    gBuffer_.sizeIs(targetBufferSize);
  }

  toc = getTime();
  log() << "-> Build buffer: " << toc - tic << " ms\n";
  tic = toc;

  // Collect data from time-slices
  HalfSliceBasisCollectorImpl::CollectedState cs = collector_->firstForwardFinalState();
  while (cs.second.vectorSize() != 0) {
    int inBufferRank = globalExchangeNumbering_->globalIndex(HalfSliceId(cs.first, HalfTimeSlice::FORWARD));
    if (inBufferRank >= 0) {
      double * targetBuffer = gBuffer_.array() + (inBufferRank * stateSize);
      bufferStateCopy(cs.second, targetBuffer); 
    }
    collector_->firstForwardFinalStateDel();
    cs = collector_->firstForwardFinalState();
  }

  cs = collector_->firstBackwardFinalState();
  std::deque<HalfSliceBasisCollectorImpl::CollectedState> localInitialStates; // Temporarily preserve initial states
  while (cs.second.vectorSize() != 0) {
    int inBufferRank = globalExchangeNumbering_->globalIndex(HalfSliceId(cs.first, HalfTimeSlice::BACKWARD));
    if (inBufferRank >= 0) {
      double * targetBuffer = gBuffer_.array() + (inBufferRank * stateSize);
      if (metric_) {
        // Pre-multiply initial local states by the metric first
        const_cast<SparseMatrix*>(metric_->stiffnessMatrix())->mult(cs.second.displacement(), targetBuffer);
        const_cast<SparseMatrix*>(metric_->massMatrix())->mult(cs.second.velocity(), targetBuffer + vectorSize_);
      } else {
        bufferStateCopy(cs.second, targetBuffer); 
      }
      localInitialStates.push_back(cs);
    }
    collector_->firstBackwardFinalStateDel();
    cs = collector_->firstBackwardFinalState();
  }

  toc = getTime();
  log() << "-> Collect data from local time-slices: " << toc - tic << " ms\n";
  tic = toc;
  
  // Setup parameters for global MPI communication
  int numCpus = mapping_->availableCpus().value(); 
  mpiParameters_.sizeIs(2 * numCpus);
  int * recv_counts = mpiParameters_.array();
  int * displacements = mpiParameters_.array() + numCpus;

  for (int i = 0; i < numCpus; ++i) {
    recv_counts[i] = globalExchangeNumbering_->stateCount(CpuRank(i)) * stateSize;
  }

  displacements[0] = 0;
  for (int i = 1; i < numCpus; ++i) {
    displacements[i] = displacements[i-1] + recv_counts[i-1];
  }

  int myCpu = localCpu_.value();
  //log() << "Allgatherv: " << recv_counts[myCpu]  << " / " << globalExchangeNumbering_->stateCount(CpuRank(myCpu)) * stateSize
  //      << " -> " << displacements[numCpus-1] + recv_counts[numCpus-1] << " / " << globalExchangeNumbering_->stateCount() * stateSize;
  timeCommunicator_->allGatherv(gBuffer_.array() + displacements[myCpu], recv_counts[myCpu], gBuffer_.array(), recv_counts, displacements);
  //log() << " complete !\n";
  
  toc = getTime();
  log() << "-> Allgather: " << toc - tic << " ms\n";
  tic = toc;

  // Add new states to projection bases
  DynamStateBasisWrapper::Ptr receivedBasis = DynamStateBasisWrapper::New(vectorSize_, globalExchangeNumbering_->stateCount(), gBuffer_.array());
  
  for (GlobalExchangeNumbering::IteratorConst it = globalExchangeNumbering_->globalIndex(); it; ++it) {
    std::pair<HalfTimeSlice::Direction, int> p = *it;
    if (p.first == HalfTimeSlice::BACKWARD) {
      // Initial (premultiplied) state
      DynamState initialOrtho(receivedBasis->state(p.second));
      metricBasis_->lastStateIs(initialOrtho);
    } else {
      // Final state
      finalBasis_->lastStateIs(receivedBasis->state(p.second));
    }
  }
  
  toc = getTime();
  log() << "-> Add new states: " << toc - tic << " ms\n";
  tic = toc;
 
  // Assemble normal matrix in parallel (local rows)
  HalfTimeSlice::Direction initStateFlag(HalfTimeSlice::BACKWARD);

  int previousMatrixSize = normalMatrix_.dim();
  int matrixSizeIncrement = globalExchangeNumbering_->stateCount(initStateFlag);
  int newMatrixSize = previousMatrixSize + matrixSizeIncrement;
  
  size_t matrixBufferSize = newMatrixSize * matrixSizeIncrement;
  if (mBuffer_.size() < matrixBufferSize) {
    mBuffer_.sizeIs(matrixBufferSize);
  }

  for (std::deque<HalfSliceBasisCollectorImpl::CollectedState>::const_iterator it = localInitialStates.begin();
      it != localInitialStates.end();
      ++it) {
    int inBufferRank = globalExchangeNumbering_->globalHalfIndex(HalfSliceId(it->first, initStateFlag));
    double * rowBuffer = mBuffer_.array() + inBufferRank * newMatrixSize;
    for (int i = 0; i < newMatrixSize; ++i) {
      rowBuffer[i] = metricBasis_->state(i) * it->second;
    }
  }
  
  toc = getTime();
  log() << "-> Assemble normal matrix: " << toc - tic << " ms\n";
  tic = toc;

  // Exchange normal matrix data  
  for (int i = 0; i < numCpus; ++i) {
    recv_counts[i] = globalExchangeNumbering_->stateCount(CpuRank(i), initStateFlag) * newMatrixSize;
  }

  displacements[0] = 0;
  for (int i = 1; i < numCpus; ++i) {
    displacements[i] = displacements[i-1] + recv_counts[i-1];
  }

  timeCommunicator_->allGatherv(mBuffer_.array() + displacements[myCpu], recv_counts[myCpu], mBuffer_.array(), recv_counts, displacements);
  
  toc = getTime();
  log() << "-> Exchange normal matrix: " << toc - tic << " ms\n";
  tic = toc;

  // Fill in normal matrix
  normalMatrix_.reSize(newMatrixSize);

  const double * originBufferBegin = mBuffer_.array();
  for (GlobalExchangeNumbering::IteratorConst it = globalExchangeNumbering_->globalHalfIndex(initStateFlag); it; ++it) {
    int rowIndex = (*it).second;
    std::copy(originBufferBegin, originBufferBegin + (rowIndex + 1), normalMatrix_[rowIndex]);
    originBufferBegin += newMatrixSize;
  }
  
  toc = getTime();
  log() << "-> Assemble normal matrix: " << toc - tic << " ms\n";
  tic = toc;

  solver_->matrixIs(normalMatrix_);
  solver_->statusIs(PivotedCholeskySolver::FACTORIZED);
   
  // HACK to reset the Reductor and Reconstructor instances 
  reductorMgr_->defaultReductionBasisIs(metricBasis_.ptr());
  reconstructorMgr_->defaultReconstructionBasisIs(finalBasis_.ptr());

  toc = getTime();
  log() << "-> Factor: " << toc - tic << " ms\n";
  tic = toc;
}

// HalfSliceCorrectionNetworkImpl Reactors

void
HalfSliceCorrectionNetworkImpl::ProjectionBuildingReactor::onStatus() {
  if (notifier()->status() == Activity::executing) {
    parent()->buildProjection();
  }
}


void
HalfSliceCorrectionNetworkImpl::SchedulingReactor::onStatus() {
  switch (notifier()->status()) {
    case Activity::free:
      // Automatic rescheduling
      notifier()->iterationIs(notifier()->currentIteration().next());
      notifier()->statusIs(Activity::scheduled);
      break;
    case Activity::scheduled:
      break;
    case Activity::executing:
      // Use the current convergence status (before it is updated) 
      parent()->globalExchangeNumbering_ = new GlobalExchangeNumbering(parent()->mapping_.ptr());
      {
        Activity * activity = parent()->projectionBuildingReactor_->notifier().ptr();
        activity->iterationIs(activity->currentIteration().next());
        activity->statusIs(Activity::scheduled);
      }
      break;
  }
}

void
HalfSliceCorrectionNetworkImpl::NonHomogeneousSchedulingReactor::onStatus() {
  if (notifier()->status() != Activity::executing ||
      notifier()->currentIteration() > IterationRank(0))
  {
    SchedulingReactor::onStatus();
  }
}
  
// HalfSliceCorrectionNetworkImpl::GlobalExchangeNumbering

HalfSliceCorrectionNetworkImpl::GlobalExchangeNumbering::GlobalExchangeNumbering(const SliceMapping * mapping)
{
  this->initialize(mapping);
}

void
HalfSliceCorrectionNetworkImpl::GlobalExchangeNumbering::initialize(const SliceMapping * mapping) {
  // Count only full timeslices
  int fullSliceCount = std::max(0, (mapping->activeSlices().value() - 1) / 2);
  int numCpus = mapping->availableCpus().value();

  // Avoid reallocations
  stateCount_.reserve(numCpus);
  initialStateCount_.reserve(numCpus);
  finalStateCount_.reserve(numCpus);
 
  stateId_.reserve(fullSliceCount * 2);
  initialStateId_.reserve(fullSliceCount);
  finalStateId_.reserve(fullSliceCount);

  // Determine location of the HalfTimeSlices updated since last correction
  typedef std::set<HalfSliceRank> HalfStateSet;
  typedef std::pair<HalfStateSet, HalfStateSet> FullStateSet; // pair<initialStates, finalStates>
  typedef std::vector<FullStateSet> TempMapping; 
  TempMapping currentStateMapping(mapping->availableCpus().value());
  
  HalfSliceRank r = HalfSliceRank(1) + mapping->convergedSlices();
  for (int fts = 0; fts < fullSliceCount; ++fts) {
    // Initial states <=> Backward
    CpuRank cpu = mapping->hostCpu(SliceId(BACKWARD_HALF_SLICE, r));
    currentStateMapping[cpu.value()].first.insert(r);
    r = r + HalfSliceCount(1); 

    // Initial states <=> Forward
    cpu = mapping->hostCpu(SliceId(FORWARD_HALF_SLICE, r));
    currentStateMapping[cpu.value()].second.insert(r);
    r = r + HalfSliceCount(1); 
  }

  // Fill the member data structures
  size_t currentFullGlobalIndex = 0;
  size_t currentInitialGlobalIndex = 0;
  size_t currentFinalGlobalIndex = 0;

  for (TempMapping::const_iterator it = currentStateMapping.begin();
       it != currentStateMapping.end();
       ++it) {
    
    initialStateCount_.push_back(it->first.size());
    finalStateCount_.push_back(it->second.size());
    stateCount_.push_back(it->first.size() + it->second.size());

    // Initial states
    for (HalfStateSet::const_iterator jt = it->first.begin(); jt != it->first.end(); ++jt) {
      HalfSliceId stateId(*jt, HalfTimeSlice::BACKWARD);
      
      initialGlobalIndex_.insert(std::make_pair(stateId, currentInitialGlobalIndex++));
      globalIndex_.insert(std::make_pair(stateId, currentFullGlobalIndex++));
      
      initialStateId_.push_back(stateId);
      stateId_.push_back(stateId);
    }
    
    // Final states
    for (HalfStateSet::const_iterator jt = it->second.begin(); jt != it->second.end(); ++jt) {
      HalfSliceId stateId(*jt, HalfTimeSlice::FORWARD);
      
      finalGlobalIndex_.insert(std::make_pair(stateId, currentFinalGlobalIndex++));
      globalIndex_.insert(std::make_pair(stateId, currentFullGlobalIndex++));
      
      finalStateId_.push_back(stateId);
      stateId_.push_back(stateId);
    }

  }

}

size_t
HalfSliceCorrectionNetworkImpl::GlobalExchangeNumbering::stateCount(HalfTimeSlice::Direction d) const {
  switch (d) {
    case HalfTimeSlice::NO_DIRECTION:
      return stateId_.size();
    case HalfTimeSlice::FORWARD:
      return finalStateId_.size();
    case HalfTimeSlice::BACKWARD:
      return initialStateId_.size();
  }
  return 0;
}

size_t
HalfSliceCorrectionNetworkImpl::GlobalExchangeNumbering::stateCount(CpuRank c, HalfTimeSlice::Direction d) const {
  try {
    switch (d) {
      case HalfTimeSlice::NO_DIRECTION:
        return stateCount_.at(c.value());
      case HalfTimeSlice::FORWARD:
        return finalStateCount_.at(c.value());
      case HalfTimeSlice::BACKWARD:
        return initialStateCount_.at(c.value());
    }
  } catch (std::out_of_range & e) {
    // Do nothing
  }
  return 0; 
}

int
HalfSliceCorrectionNetworkImpl::GlobalExchangeNumbering::globalIndex(const HalfSliceId & id) const {
  IndexMap::const_iterator it = globalIndex_.find(id);
  return (it != globalIndex_.end()) ? static_cast<int>(it->second) : -1; 
}

HalfSliceCorrectionNetworkImpl::GlobalExchangeNumbering::IteratorConst
HalfSliceCorrectionNetworkImpl::GlobalExchangeNumbering::globalIndex() const {
  return IteratorConst(this->globalIndex_.begin(), this->globalIndex_.end());
}

int
HalfSliceCorrectionNetworkImpl::GlobalExchangeNumbering::globalHalfIndex(const HalfSliceId & id) const {
  int result = -1;
  
  IndexMap::const_iterator it;
  switch (id.direction()) {
    case HalfTimeSlice::NO_DIRECTION:
      break;
    case HalfTimeSlice::FORWARD:
      it = finalGlobalIndex_.find(id);
      if (it != finalGlobalIndex_.end())
        result = static_cast<int>(it->second);
      break;
    case HalfTimeSlice::BACKWARD:
      it = initialGlobalIndex_.find(id);
      if (it != initialGlobalIndex_.end())
        result = static_cast<int>(it->second);
  }

  return result;
}

HalfSliceCorrectionNetworkImpl::GlobalExchangeNumbering::IteratorConst
HalfSliceCorrectionNetworkImpl::GlobalExchangeNumbering::globalHalfIndex(HalfTimeSlice::Direction d) const {
  switch (d) {
    case HalfTimeSlice::NO_DIRECTION:
      break;
    case HalfTimeSlice::FORWARD:
      return IteratorConst(this->finalGlobalIndex_.begin(), this->finalGlobalIndex_.end());
    case HalfTimeSlice::BACKWARD:
      return IteratorConst(this->initialGlobalIndex_.begin(), this->initialGlobalIndex_.end());
  }
  return IteratorConst();
}

HalfSliceId
HalfSliceCorrectionNetworkImpl::GlobalExchangeNumbering::stateId(int gfi) const {
  try {
    return stateId_.at(gfi);
  } catch (std::out_of_range & e) {
    // Do nothing
  }
  return HalfSliceId();
}

HalfSliceId
HalfSliceCorrectionNetworkImpl::GlobalExchangeNumbering::stateId(int ghi, HalfTimeSlice::Direction d) const {
  try {
    switch (d) {
      case HalfTimeSlice::NO_DIRECTION:
        break;
      case HalfTimeSlice::FORWARD:
        return finalStateId_.at(ghi);
      case HalfTimeSlice::BACKWARD:
        return initialStateId_.at(ghi);
    }
  } catch (std::out_of_range & e) {
    // Do nothing
  }
  return HalfSliceId();
}

} // end namespace Pita
