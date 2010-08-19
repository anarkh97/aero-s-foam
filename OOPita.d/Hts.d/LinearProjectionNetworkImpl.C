#include "LinearProjectionNetworkImpl.h"

#include "../DynamStateBasisWrapper.h"

#include "../DynamStateOps.h"
#include <Comm.d/Communicator.h>

#include <stdexcept>

#include <Timers.d/GetTime.h>

namespace Pita { namespace Hts {

LinearProjectionNetworkImpl::LinearProjectionNetworkImpl(size_t vSize,
                                             Communicator * timeComm,
                                             CpuRank myCpu,
                                             const SliceMapping * mapping,
                                             BasisCollectorImpl * collector,
                                             const DynamOps * metric,
                                             RankDeficientSolver * solver) :
  vectorSize_(vSize),
  timeCommunicator_(timeComm),
  localCpu_(myCpu),
  mapping_(mapping),
  metric_(metric),
  gBuffer_(),
  mBuffer_(),
  mpiParameters_(),
  localBasis_(),
  metricBasis_(DynamStatePlainBasis::New(vectorSize_)),
  finalBasis_(DynamStatePlainBasis::New(vectorSize_)),
  originalMetricBasis_(DynamStatePlainBasis::New(vectorSize_)),
  originalFinalBasis_(DynamStatePlainBasis::New(vectorSize_)),
  normalMatrix_(),
  transmissionMatrix_(),
  reprojectionMatrix_(),
  solver_(solver),
  collector_(collector),
  globalExchangeNumbering_()
{}

void
LinearProjectionNetworkImpl::prepareProjection() {
  GlobalExchangeNumbering::Ptr numbering = new GlobalExchangeNumbering(mapping_.ptr());
  globalExchangeNumbering_.push_back(numbering);
}

void
LinearProjectionNetworkImpl::buildProjection() {
#ifndef NDEBUG
  double tic = getTime(), toc;
#endif /* NDEBUG*/

  int previousMatrixSize = normalMatrix_.dim();

  // Build buffer
  size_t stateSize = 2 * vectorSize_;
  size_t targetBufferSize = stateSize * globalExchangeNumbering_.back()->stateCount();
  if (gBuffer_.size() < targetBufferSize) {
    gBuffer_.sizeIs(targetBufferSize);
  }

#ifndef NDEBUG
  toc = getTime();
  log() << "      -> Build buffer: " << toc - tic << " ms\n";
  tic = toc;
#endif /* NDEBUG*/
  
  // Collect data from time-slices
  // 1) Final states
  BasisCollectorImpl::CollectedState cs = collector_->firstForwardFinalState();
  while (cs.second.vectorSize() != 0) {
    int inBufferRank = globalExchangeNumbering_.back()->globalIndex(HalfSliceId(cs.first, HalfTimeSlice::FORWARD));
    if (inBufferRank >= 0) {
      double * targetBuffer = gBuffer_.array() + (inBufferRank * stateSize);
      bufferStateCopy(cs.second, targetBuffer); 
      //log() << "Final state # " << inBufferRank << " disp[2] = " << targetBuffer[2] << "\n";
      int accumulatedIndex = inBufferRank + (2 * previousMatrixSize);
      localBasis_.insert(std::make_pair(accumulatedIndex, cs.second));
    }
    collector_->firstForwardFinalStateDel();
    cs = collector_->firstForwardFinalState();
  }

  // 2) Initial states
  cs = collector_->firstBackwardFinalState();
  while (cs.second.vectorSize() != 0) {
    int inBufferRank = globalExchangeNumbering_.back()->globalIndex(HalfSliceId(cs.first, HalfTimeSlice::BACKWARD));
    if (inBufferRank >= 0) {
      double * targetBuffer = gBuffer_.array() + (inBufferRank * stateSize);
      if (metric_) {
        // Pre-multiply initial local states by the metric first
        mult(metric_.ptr(), cs.second, targetBuffer);
      } else {
        log() << "Warning, no metric found !\n";
        bufferStateCopy(cs.second, targetBuffer); 
      }
      //log() << "Metric state # " << inBufferRank << " disp[2] = " << targetBuffer[2] << "\n";
      int accumulatedIndex = inBufferRank + (2 * previousMatrixSize);
      localBasis_.insert(std::make_pair(accumulatedIndex, cs.second));
    }
    collector_->firstBackwardFinalStateDel();
    cs = collector_->firstBackwardFinalState();
  }

#ifndef NDEBUG
  toc = getTime();
  log() << "      -> Collect data from local time-slices: " << toc - tic << " ms\n";
  tic = toc;
#endif /* NDEBUG*/
  
  // Setup parameters for global MPI communication
  int numCpus = mapping_->availableCpus().value(); 
  mpiParameters_.sizeIs(2 * numCpus);
  int * recv_counts = mpiParameters_.array();
  int * displacements = mpiParameters_.array() + numCpus;

  for (int i = 0; i < numCpus; ++i) {
    recv_counts[i] = globalExchangeNumbering_.back()->stateCount(CpuRank(i)) * stateSize;
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
  
#ifndef NDEBUG
  toc = getTime();
  log() << "      -> Allgather: " << toc - tic << " ms\n";
  tic = toc;
#endif /* NDEBUG*/

  // Add new states to projection bases
  DynamStateBasisWrapper::Ptr receivedBasis = DynamStateBasisWrapper::New(vectorSize_, globalExchangeNumbering_.back()->stateCount(), gBuffer_.array());

  /*for (int i = 0; i < globalExchangeNumbering_.back()->stateCount(); ++i) {
    log() << "Received State # " << i << " disp[2] = " << receivedBasis->state(i).displacement()[2] << "\n";
  }*/

  for (GlobalExchangeNumbering::IteratorConst it = globalExchangeNumbering_.back()->globalIndex(); it; ++it) {
    std::pair<HalfTimeSlice::Direction, int> p = *it;
    DynamStatePlainBasis * targetBasis = (p.first == HalfTimeSlice::BACKWARD) ?
                                        originalMetricBasis_.ptr() : // Initial state (premultiplied)
                                        originalFinalBasis_.ptr();   // Final state
    targetBasis->lastStateIs(receivedBasis->state(p.second));
  }
 
  /*for (int i = 0; i < finalBasis_->stateCount(); ++i) {
    log() << "Final State # " << i << " disp[2] = " << finalBasis_->state(i).displacement()[2] << "\n";
  }

  for (int i = 0; i < metricBasis_->stateCount(); ++i) {
    log() << "Metric State # " << i << " disp[2] = " << metricBasis_->state(i).displacement()[2] << "\n";
  }*/

#ifndef NDEBUG
  toc = getTime();
  log() << "      -> Add new states: " << toc - tic << " ms\n";
  tic = toc;
#endif /* NDEBUG*/
 
  // Assemble normal and reprojection matrices in parallel (local rows)
  int matrixSizeIncrement = globalExchangeNumbering_.back()->stateCount(HalfTimeSlice::BACKWARD);
  int newMatrixSize = previousMatrixSize + matrixSizeIncrement;
 
  log() << "*** Matrix size: previous = " << previousMatrixSize << ", increment = " << matrixSizeIncrement << ", newSize = " << newMatrixSize << "\n";
  
  size_t matrixBufferSize = (2 * newMatrixSize) * newMatrixSize; // For normalMatrix and reprojectionMatrix data
  if (mBuffer_.size() < matrixBufferSize) {
    mBuffer_.sizeIs(matrixBufferSize);
  }
  for (int i = 0; i < numCpus; ++i) {
    recv_counts[i] = 0;
    for (NumberingList::const_iterator it = globalExchangeNumbering_.begin(); it != globalExchangeNumbering_.end(); ++it) {
      recv_counts[i] += (*it)->stateCount(CpuRank(i));
    }
    recv_counts[i] *= newMatrixSize;
  }

  displacements[0] = 0;
  for (int i = 1; i < numCpus; ++i) {
    displacements[i] = displacements[i-1] + recv_counts[i-1];
  }
  
  // TODO Better efficiency: Do not recompute every coefs
  double * rowBuffer = mBuffer_.array() + displacements[myCpu];
  for (LocalBasis::const_iterator it = localBasis_.begin();
      it != localBasis_.end(); ++it) {
    for (int i = 0; i < newMatrixSize; ++i) {
      rowBuffer[i] = originalMetricBasis_->state(i) * it->second;
    }
    rowBuffer += newMatrixSize;
  }
  
#ifndef NDEBUG
  toc = getTime();
  log() << "      -> Compute normal and reprojection matrices: " << toc - tic << " ms\n";
  tic = toc;
#endif /* NDEBUG*/

  // Exchange normal/reprojection matrix data  
  timeCommunicator_->allGatherv(mBuffer_.array() + displacements[myCpu], recv_counts[myCpu], mBuffer_.array(), recv_counts, displacements);
  
#ifndef NDEBUG
  toc = getTime();
  log() << "      -> Exchange normal matrix: " << toc - tic << " ms\n";
  tic = toc;
#endif /* NDEBUG*/

  /*log() << "Buffer is\n";
  for (int i = 0; i < newMatrixSize; ++i) {
    for (int j = 0; j < 2 * newMatrixSize; ++j) {
      log() << "[" << i << "," << j << "] " << mBuffer_[i * newMatrixSize + j] << "\n";
    }
  }*/

  // Assemble updated normal & reprojection matrices
  normalMatrix_.reSize(newMatrixSize);
  transmissionMatrix_.reSize(newMatrixSize);

  // Loop on iterations
  int originRowIndex = 0; // Target matrix base row index for current iteration 
  for (NumberingList::const_iterator it = globalExchangeNumbering_.begin(); it != globalExchangeNumbering_.end(); ++it) {
    GlobalExchangeNumbering::PtrConst numbering = *it;
    GlobalExchangeNumbering::IteratorConst jt_i = numbering->globalHalfIndex(HalfTimeSlice::BACKWARD);
    GlobalExchangeNumbering::IteratorConst jt_f = numbering->globalHalfIndex(HalfTimeSlice::FORWARD);
   
    //log() << "[States from iteration # " << std::distance((NumberingList::const_iterator)(globalExchangeNumbering_.begin()), it) << "]\n";

    // Loop on cpus
    for (int cpu = 0; cpu < numCpus; ++cpu) {
      //log() << "Cpu # " << cpu << "\n";
      double * originBufferBegin = mBuffer_.array() + displacements[cpu]; // Position in AllgatherBuffer
      
      int initialStateCountInIter = numbering->stateCount(CpuRank(cpu), HalfTimeSlice::BACKWARD);
      //log() << "# initial states on cpu = " << initialStateCountInIter << "\n";
      for (int s = 0; s < initialStateCountInIter; ++s) {
        int rowIndex = (*jt_i).second + originRowIndex; // Row in target normal matrix
        //log() << "Normal matrix row # = " << rowIndex << "\n";
        //log() << "Buffer row # = " << std::distance(mBuffer_.array(), originBufferBegin) / newMatrixSize << "\n";
        //log() << "State id = " << numbering->stateId(std::distance(mBuffer_.array(), originBufferBegin) / newMatrixSize) << "\n";
        std::copy(originBufferBegin, originBufferBegin + newMatrixSize, normalMatrix_[rowIndex]);
        originBufferBegin += newMatrixSize;
        ++jt_i;
      }

      int finalStateCountInIter = numbering->stateCount(CpuRank(cpu), HalfTimeSlice::FORWARD);
      //log() << "# final states on cpu = " << finalStateCountInIter << "\n";
      for (int s = 0; s < finalStateCountInIter; ++s) {
        int rowIndex = (*jt_f).second + originRowIndex; // Row in target reprojection matrix
        //log() << "Reprojection matrix row # = " << rowIndex << "\n";
        //log() << "Buffer row # = " << std::distance(mBuffer_.array(), originBufferBegin) / newMatrixSize << "\n";
        //log() << "State id = " << numbering->stateId(std::distance(mBuffer_.array(), originBufferBegin) / newMatrixSize) << "\n";
        std::copy(originBufferBegin, originBufferBegin + newMatrixSize, transmissionMatrix_[rowIndex]);
        originBufferBegin += newMatrixSize;
        ++jt_f;
      }
      // Update buffer cpu base displacement for next iteration
      displacements[cpu] += newMatrixSize * (initialStateCountInIter + finalStateCountInIter);
    } 
    
    originRowIndex += numbering->stateCount(HalfTimeSlice::BACKWARD); // Udpate target matrix base row for next iteration
  }

  /*log() << "ReprojectionMatrix:\n";
  for (int i = 0; i < newMatrixSize; ++i) {
    for (int j = 0; j < newMatrixSize; ++j) {
      log() << reprojectionMatrix_[i][j] << " ";
    }
    log() << "\n";
  } 
  
  log() << "normalMatrix:\n";
  for (int i = 0; i < newMatrixSize; ++i) {
    for (int j = 0; j < newMatrixSize; ++j) {
      log() << normalMatrix_[i][j] << " ";
    }
    log() << "\n";
  } */

#ifndef NDEBUG
  toc = getTime();
  log() << "      -> Assemble normal matrix: " << toc - tic << " ms\n";
  tic = toc;
#endif /* NDEBUG*/

  //log() << "Check symmetry\n";
  
  // Check x^T (M y) = y^T (M x) and x^T (K y) = y^T (M x)
  
  /*DynamState x = finalBasis_->state(0);
  DynamState y = finalBasis_->state(1); 

  log() << x.vectorSize() << "\n"; 
  
  // M
  Vector x_v = x.velocity();
  log() << x_v.size() << "\n"; 
  Vector y_v = y.velocity();
  log() << y_v.size() << "\n"; 
  Vector My_v(y_v.size());
  log() << My_v.size() << "\n";
  
  const_cast<SparseMatrix*>(metric_->massMatrix())->mult(y_v, My_v); 
  double x_vTMy_v = x_v * My_v;
  Vector Mx_v(x_v.size());
  log() << Mx_v.size() << "\n";
  const_cast<SparseMatrix*>(metric_->massMatrix())->mult(x_v, Mx_v); 
  double y_vTMx_v = y_v * Mx_v;

  log() << " x^T (M y) / y^T (M x) = " << x_vTMy_v << " / " << y_vTMx_v << "\n";
  
  // K
  Vector x_d = x.displacement();
  Vector y_d = y.displacement();
  
  Vector Ky_d(y_d.size());
  const_cast<SparseMatrix*>(metric_->stiffnessMatrix())->mult(y_d, Ky_d); 
  double x_dTKy_d = x_d * Ky_d;
  Vector Kx_d(x_d.size());
  const_cast<SparseMatrix*>(metric_->stiffnessMatrix())->mult(x_d, Kx_d); 
  double y_dTKx_d = y_d * Kx_d;

  log() << " x^T (K y) / y^T (K x) = " << x_dTKy_d << " / " << y_dTKx_d << "\n";*/

  solver_->transposedMatrixIs(normalMatrix_);

  /*reprojectionMatrix_.copy(transmissionMatrix_);
  metricBasis_->stateBasisDel();
  finalBasis_->stateBasisDel();
  metricBasis_->lastStateBasisIs(originalMetricBasis_.ptr());
  finalBasis_->lastStateBasisIs(originalFinalBasis_.ptr());*/

  reprojectionMatrix_.reSize(solver_->factorRank());
  metricBasis_->stateBasisDel();
  finalBasis_->stateBasisDel();

  for (int compactIndex = 0; compactIndex < solver_->factorRank(); ++compactIndex) {
    int originalIndex = solver_->factorPermutation(compactIndex);

    for (int j = 0; j < solver_->factorRank(); ++j) {
      reprojectionMatrix_[compactIndex][j] = transmissionMatrix_[originalIndex][solver_->factorPermutation(j)];
    }

    metricBasis_->lastStateIs(originalMetricBasis_->state(originalIndex));
    finalBasis_->lastStateIs(originalFinalBasis_->state(originalIndex));
  }

  /*log() << "TransmissionMatrix (size = " << transmissionMatrix_.dim() << "):\n";
  for (int i = 0; i < transmissionMatrix_.dim(); ++i) {
    for (int j = 0; j < transmissionMatrix_.dim(); ++j) {
      log() << transmissionMatrix_[i][j] << " ";
    }
    log() << "\n";
  } 

  log() << "ReprojectionMatrix (size = " << reprojectionMatrix_.dim() << "):\n";
  for (int i = 0; i < reprojectionMatrix_.dim(); ++i) {
    for (int j = 0; j < reprojectionMatrix_.dim(); ++j) {
      log() << reprojectionMatrix_[i][j] << " ";
    }
    log() << "\n";
  }
  
  log() << "Permutation =";
  for (int i = 0; i < solver_->factorRank(); ++i) {
    log() << " " << solver_->factorPermutation(i);
  }
  log() << "\n";

  log() << "check\n";
  for (int i = 0; i < reprojectionMatrix_.dim(); ++i) {
    for (int j = 0; j < reprojectionMatrix_.dim(); ++j) {
      log() << reprojectionMatrix_[i][j] - transmissionMatrix_[solver_->factorPermutation(i)][solver_->factorPermutation(j)] << " ";
    }
    log() << "\n";
  }*/
  
  solver_->orderingIs(RankDeficientSolver::COMPACT);

  /*log() << "Compact permutation =";
  for (int i = 0; i < solver_->factorRank(); ++i) {
    log() << " " << solver_->factorPermutation(i);
  }
  log() << "\n";*/

  /*if (timeCommunicator_->myID() == 0) {
    log() << "Matrix\n";
    for (int i = 0; i < solver_->matrixSize(); ++i) {
      for (int j = 0; j < solver_->matrixSize(); ++j) {
        log() << const_cast<FullSquareMatrix &>(solver_->transposedMatrix())[i][j] << " ";
      }
      log() << "\n";
    }
  }*/
  
  //solver_->statusIs(RankDeficientSolver::FACTORIZED);
  
  /*if (timeCommunicator_->myID() == 0) {
    log() << "Factor\n";
    for (int i = 0; i < solver_->factorRank(); ++i) {
      for (int j = 0; j < solver_->factorRank(); ++j) {
        log() << const_cast<FullSquareMatrix &>(solver_->transposedMatrix())[i][j] << " ";
      }
      log() << "\n";
    }
  }*/
  
#ifndef NDEBUG
  toc = getTime();
  log() << "      -> Factor: " << toc - tic << " ms\n";
  tic = toc;
#endif /* NDEBUG*/

  // Debug
  //log() << "Debug\n";
  log() << "*** Projector rank = " << solver_->factorRank() << "\n";
  /*log() << "Permutation =";
  for (int i = 0; i < solver_->factorRank(); ++i) {
    log() << " " << solver_->factorPermutation(i);
  }
  log() << "\n";*/
  
}
  
// LinearProjectionNetworkImpl::GlobalExchangeNumbering

LinearProjectionNetworkImpl::GlobalExchangeNumbering::GlobalExchangeNumbering(const SliceMapping * mapping)
{
  this->initialize(mapping);
}

void
LinearProjectionNetworkImpl::GlobalExchangeNumbering::initialize(const SliceMapping * mapping) {
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
    CpuRank cpu = mapping->hostCpu(r);
    currentStateMapping[cpu.value()].first.insert(r);
    r = r + HalfSliceCount(1); 

    // Initial states <=> Forward
    cpu = mapping->hostCpu(r);
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
LinearProjectionNetworkImpl::GlobalExchangeNumbering::stateCount() const {
  return stateId_.size();
}

size_t
LinearProjectionNetworkImpl::GlobalExchangeNumbering::stateCount(HalfTimeSlice::Direction d) const {
  switch (d) {
    case HalfTimeSlice::NO_DIRECTION:
      return 0;
    case HalfTimeSlice::FORWARD:
      return finalStateId_.size();
    case HalfTimeSlice::BACKWARD:
      return initialStateId_.size();
  }
  throw Fwk::InternalException("in LinearProjectionNetworkImpl::GlobalExchangeNumbering::stateCount");
}

size_t
LinearProjectionNetworkImpl::GlobalExchangeNumbering::stateCount(CpuRank c) const {
  try {
    return stateCount_.at(c.value());
  } catch (std::out_of_range & e) {
    // Do nothing
  }
  return 0; 
}

size_t
LinearProjectionNetworkImpl::GlobalExchangeNumbering::stateCount(CpuRank c, HalfTimeSlice::Direction d) const {
  try {
    switch (d) {
      case HalfTimeSlice::NO_DIRECTION:
        return 0;
      case HalfTimeSlice::FORWARD:
        return finalStateCount_.at(c.value());
      case HalfTimeSlice::BACKWARD:
        return initialStateCount_.at(c.value());
      default:
        throw Fwk::InternalException("in LinearProjectionNetworkImpl::GlobalExchangeNumbering::stateCount");
    }
  } catch (std::out_of_range & e) {
    // Do nothing
  }
  return 0;
}

int
LinearProjectionNetworkImpl::GlobalExchangeNumbering::globalIndex(const HalfSliceId & id) const {
  IndexMap::const_iterator it = globalIndex_.find(id);
  return (it != globalIndex_.end()) ? static_cast<int>(it->second) : -1; 
}

LinearProjectionNetworkImpl::GlobalExchangeNumbering::IteratorConst
LinearProjectionNetworkImpl::GlobalExchangeNumbering::globalIndex() const {
  return IteratorConst(this->globalIndex_.begin(), this->globalIndex_.end());
}

int
LinearProjectionNetworkImpl::GlobalExchangeNumbering::globalHalfIndex(const HalfSliceId & id) const {
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

LinearProjectionNetworkImpl::GlobalExchangeNumbering::IteratorConst
LinearProjectionNetworkImpl::GlobalExchangeNumbering::globalHalfIndex(HalfTimeSlice::Direction d) const {
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
LinearProjectionNetworkImpl::GlobalExchangeNumbering::stateId(int gfi) const {
  try {
    return stateId_.at(gfi);
  } catch (std::out_of_range & e) {
    // Do nothing
  }
  return HalfSliceId();
}

HalfSliceId
LinearProjectionNetworkImpl::GlobalExchangeNumbering::stateId(int ghi, HalfTimeSlice::Direction d) const {
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

} /* end namespace Hts */ } /* end namespace Pita */
