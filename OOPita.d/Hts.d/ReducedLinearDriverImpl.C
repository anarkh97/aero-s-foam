#include "Fwk.h"
#include "Types.h"

#include "ReducedLinearDriverImpl.h"

#include "../DynamState.h"
#include <Problems.d/DynamDescr.h>

#include "../LinearGenAlphaIntegrator.h"
#include "../HomogeneousGenAlphaIntegrator.h"
#include "../IntegratorPropagator.h"

#include "SliceMapping.h"

#include "../RankDeficientSolver.h"
#include "../NearSymmetricSolver.h"
#include "../PivotedCholeskySolver.h"

#include "BasisCollectorImpl.h"
#include "NonHomogeneousBasisCollectorImpl.h"
#include "CorrectionNetworkImpl.h"

#include "../Seed.h"

#include "../RemoteStateMpiImpl.h"

#include "../IncrementalPostProcessor.h"
#include "LinearFineIntegratorManager.h"
#include "PropagatorManager.h"
#include "HalfTimeSliceImpl.h"

#include "ReducedFullTimeSliceImpl.h"
#include "../JumpProjectorImpl.h"
#include "../UpdatedSeedAssemblerImpl.h"
#include "LocalCorrectionTimeSlice.h"

#include "../IntegratorSeedInitializer.h"
#include "RemoteSeedInitializerServer.h"
#include "../RemoteSeedInitializerProxy.h"

#include "../RemoteDynamPropagatorProxy.h"
#include "RemoteCoarseCorrectionServer.h"

#include "../CommSplitter.h"
#include <Comm.d/Communicator.h>

#include <Driver.d/Domain.h>

#include "LocalNetwork.h"

namespace Pita { namespace Hts {

ReducedLinearDriverImpl::ReducedLinearDriverImpl(SingleDomainDynamic * pbDesc,
                                                 Domain * domain,
                                                 SolverInfo * solverInfo,
                                                 Communicator * baseComm) :
  probDesc_(pbDesc),
  domain_(domain),
  solverInfo_(solverInfo),
  baseComm_(baseComm)
{}

// Main routine
void
ReducedLinearDriverImpl::solve() {
  log() << "Begin Reversible Linear Pita\n";
  double tic = getTime(); // Total time

  preprocess();

  /* Summarize problem and parameters */
  log() << "\n"; 
  log() << "Slices = " << mapping_->totalSlices() << ", MaxActive = " << mapping_->maxWorkload() << ", Cpus = " << mapping_->availableCpus() << "\n";
  log() << "Iteration count = " << lastIteration_ << "\n"; 
  log() << "dt = " << fineTimeStep_ << ", J/2 = " << halfSliceRatio_ << ", Dt = J*dt = " << coarseTimeStep_ << ", Tf = Slices*(J/2)*dt = " << finalTime_ << "\n";
  if (noForce_) { log() << "No external force\n"; }
  if (remoteCoarse_) { log() << "Remote coarse time-integration\n"; }
  log() << "VectorSize = " << vectorSize_ << " dofs\n";
  log() << "Projector tol = " << projectorTolerance_ << "\n";

  /* Determine process task */
  if (remoteCoarse_) {
    CommSplitter::Ptr commSplitter = CommSplitter::New(baseComm(), CpuCount(1)); // Specializes one cpu

    if (commSplitter->localGroup() == CommSplitter::STANDARD) {
      solveParallel(commSplitter->splitComm(), commSplitter->interComm());
    } else {
      solveCoarse(commSplitter->interComm());
    }
  } else {
    solveParallel(baseComm(), NULL);
  }
  
  double toc = getTime();
  log() << "\n" << "Total Time = " << (toc - tic) / 1000.0 << " s\n";
  log() << "\n" << "End Reversible Linear Pita\n";
}

void
ReducedLinearDriverImpl::preprocess() {
  double tic = getTime();
  
  probDesc()->preProcess();
  
  /* Space-domain */
  vectorSize_ = probDesc()->solVecInfo();
  dynamOpsMgr_= LinearDynamOps::Manager::New(probDesc());

  /* Time-domain */
  fineTimeStep_ = Seconds(solverInfo()->getTimeStep());
  halfSliceRatio_ = TimeStepCount(solverInfo()->Jratio / 2);
  sliceRatio_ = TimeStepCount(halfSliceRatio_.value() * 2);
  coarseTimeStep_ = fineTimeStep_ * sliceRatio_.value(); 
  initialTime_ = Seconds(solverInfo()->initialTime);
  finalTime_ = Seconds(solverInfo()->tmax);
 
  HalfSliceCount numSlices(static_cast<int>(ceil((finalTime_.value() - initialTime_.value()) / (halfSliceRatio_.value() * fineTimeStep_.value()))));
  FullSliceCount fullTimeSlices((numSlices.value() / 2) + (numSlices.value() % 2));
  numSlices = HalfSliceCount(fullTimeSlices.value() * 2);
  finalTime_ = fineTimeStep_ * Seconds(numSlices.value() * halfSliceRatio_.value()); // To have a whole number of full time-slices 
  
  /* Main options */
  noForce_ = solverInfo()->NoForcePita;
  remoteCoarse_ = solverInfo()->remoteCoarse && (baseComm()->numCPUs() > 1);

  /* Load balancing */ 
  CpuCount numCpus(baseComm()->numCPUs() - (remoteCoarse_ ? 1 : 0));
  HalfSliceCount maxActive(solverInfo()->numTSperCycleperCPU);
  mapping_ = SliceMapping::New(fullTimeSlices, numCpus, maxActive);

  /* Other parameters */ 
  lastIteration_ = IterationRank(solverInfo()->kiter);
  projectorTolerance_ = solverInfo()->pitaProjTol;
  coarseRhoInfinity_ = 1.0; // TODO Could be set in input file
 
  double toc = getTime();
  log() << "\n";
  log() << "Total Preprocessing Time = " << (toc - tic) / 1000.0 << " s\n";
}


void
ReducedLinearDriverImpl::solveParallel(Communicator * timeComm, Communicator * coarseComm) {
  log() << "\n";
  log() << "Parallel solver initialization\n";

  /* Local process identity */
  CpuRank myCpu(timeComm->myID());

  /* Fine time integration operators */
  const double noDampingRhoInfinity = 1.0; // No damping allowed for time-reversible
  GeneralizedAlphaParameter fineIntegrationParam(fineTimeStep_, noDampingRhoInfinity);
  LinearDynamOps::Ptr dynamOps = dynamOpsMgr_->dynOpsNew(fineIntegrationParam); // Expensive operation 
  
  /* Post processing */ 
  PostProcessing::Manager::Ptr postProcessingMgr = buildPostProcessor(myCpu);

  /* Local Basis Collector */
  BasisCollectorImpl::Ptr collector = buildBasisCollector(); 

  /* Correction */
  RankDeficientSolver::Ptr normalMatrixSolver = NearSymmetricSolver::New(projectorTolerance_); // TODO Other solver implementations ?
  
  CorrectionNetworkImpl::Ptr correctionMgr = CorrectionNetworkImpl::New(
      vectorSize_,
      timeComm,
      myCpu,
      mapping_.ptr(),
      collector.ptr(),
      dynamOps.ptr(),
      normalMatrixSolver.ptr());

  /* HalfTimeSlices and other local tasks */
  HalfTimeSlice::Manager::Ptr hsMgr = buildHalfTimeSliceManager(fineIntegrationParam, postProcessingMgr.ptr(), collector.ptr());

  JumpProjector::Manager::Ptr jpMgr = JumpProjectorImpl::Manager::New(correctionMgr->projectionBasis());
  ReducedFullTimeSlice::Manager::Ptr fsMgr = ReducedFullTimeSliceImpl::Manager::New(correctionMgr->reprojectionMatrix(), correctionMgr->normalMatrixSolver());
  UpdatedSeedAssembler::Manager::Ptr usaMgr = UpdatedSeedAssemblerImpl::Manager::New(correctionMgr->propagatedBasis());

  RemoteState::MpiManager::Ptr commMgr = RemoteState::MpiManager::New(timeComm);
  commMgr->vectorSizeIs(vectorSize_);

  /* Coarse time-integration */ 
  CorrectionTimeSlice::Manager::Ptr ctsMgr = buildCoarseCorrection(coarseComm);

  /* Local tasks */
  LocalNetwork::Ptr network = new LocalNetwork(
      mapping_.ptr(),
      myCpu,
      hsMgr.ptr(),
      jpMgr.ptr(),
      fsMgr.ptr(),
      usaMgr.ptr(),
      ctsMgr.ptr(),
      commMgr.ptr());

  /* Initial seeds */
  IterationRank iteration;
 
  if (noForce_) {
    SeedInitializer::Ptr seedInitializer = buildSeedInitializer(!remoteCoarse_, coarseComm); 
    iteration = IterationRank(0);
    
    LocalNetwork::MainSeedMap mainSeeds = network->activeMainSeeds();
    for (LocalNetwork::MainSeedMap::iterator it = mainSeeds.begin();
        it != mainSeeds.end();
        ++it) {
      log() << "Initial Seed " << it->first << " " << it->second->name() << "\n";
      DynamState initSeed = seedInitializer->initialSeed(it->first);
      Seed::Ptr seed = it->second;
      seed->stateIs(initSeed);
      seed->iterationIs(iteration);
      seed->statusIs(it->first == SliceRank(0) ? Seed::CONVERGED : Seed::ACTIVE);
    }

    /* Begin iterations */
    log() << "\nIteration " << iteration << "\n";

  } else {
    iteration = IterationRank(-1);

    LocalNetwork::SeedMap mainSeeds = network->mainSeeds();
    for (LocalNetwork::SeedMap::iterator it = mainSeeds.begin();
        it != mainSeeds.end();
        ++it) {
      log() << "Zero initial Seed " << it->first << " " << it->second->name() << "\n";
      Seed::Ptr seed = it->second;
      seed->stateIs(seed->name() == "M0" ? initialSeed() : DynamState(vectorSize_, 0.0));
      seed->iterationIs(iteration);
      Seed::Status status = (it->first == HalfSliceRank(0)) ? Seed::CONVERGED : Seed::ACTIVE;
      if (it->first.value() % 2 != 0) {
        status = Seed::SPECIAL;
      }
      seed->statusIs(status);
    }
   
    log() << "\nAffine term precomputation\n";
    LocalNetwork::TaskList halfSlices = network->halfTimeSlices();
    for (LocalNetwork::TaskList::iterator it = halfSlices.begin();
        it != halfSlices.end();
        ++it) {
      NamedTask::Ptr task = *it;
      log() << "Task: " << task->name() << "\n";
      task->iterationIs(iteration);
    }

    if (iteration < lastIteration_) {

      network->convergedSlicesInc();

      iteration = iteration.next(); 
      log() << "\nIteration " << iteration << "\n";
      // Propagated Seed Synchronization 
      LocalNetwork::TaskList leftSeedSyncs = network->activeLeftSeedSyncs();
      for (LocalNetwork::TaskList::iterator it = leftSeedSyncs.begin();
          it != leftSeedSyncs.end();
          ++it) {
        NamedTask::Ptr task = *it;
        log() << "Task: " << task->name() << "\n";
        task->iterationIs(iteration);
      }

      // Coarse propagation
      LocalNetwork::TaskList coarseCorr = network->activeCoarseTimeSlices();
      for (LocalNetwork::TaskList::iterator it = coarseCorr.begin();
          it != coarseCorr.end();
          ++it) {
        NamedTask::Ptr task = *it;
        log() << "Task: " << task->name() << "\n";
        task->iterationIs(iteration);
      }

      // Correction Synchronization 
      LocalNetwork::TaskList synchronizationTasks = network->activeFullCorrectionSyncs();
      for (LocalNetwork::TaskList::iterator it = synchronizationTasks.begin();
          it != synchronizationTasks.end();
          ++it) {
        NamedTask::Ptr task = *it;
        log() << "Task: " << task->name() << "\n";
        task->iterationIs(iteration);
      }

      // Update Seeds
      LocalNetwork::TaskList assemblyTasks = network->activeSeedAssemblers();
      for (LocalNetwork::TaskList::iterator it = assemblyTasks.begin();
          it != assemblyTasks.end();
          ++it) {
        NamedTask::Ptr task = *it;
        log() << "Task: " << task->name() << "\n";
        task->iterationIs(iteration);
      }
    }
  }
  
  // Local Fine Propagation
  if (iteration < lastIteration_) {
    LocalNetwork::TaskList halfSlices = network->activeHalfTimeSlices();
    for (LocalNetwork::TaskList::iterator it = halfSlices.begin();
        it != halfSlices.end();
        ++it) {
      NamedTask::Ptr task = *it;
      log() << "Task: " << task->name() << "\n";
      task->iterationIs(iteration);
    }
  }

  // Iteration loop
  for (iteration = IterationRank(1); iteration < lastIteration_; iteration = iteration.next()) {

    log() << "\nIteration " << iteration << "\n";

    // Projection Basis
    correctionMgr->buildProjection();
    size_t reducedBasisSize = correctionMgr->reducedBasisSize();
    log() << "ReducedBasisSize = " << reducedBasisSize << "\n";
    commMgr->reducedStateSizeIs(reducedBasisSize);

    // Next iteration 
    network->convergedSlicesInc();

    // Propagated Seed Synchronization 
    LocalNetwork::TaskList leftSeedSyncs = network->activeLeftSeedSyncs();
    for (LocalNetwork::TaskList::iterator it = leftSeedSyncs.begin();
        it != leftSeedSyncs.end();
        ++it) {
      NamedTask::Ptr task = *it;
      log() << "Task: " << task->name() << "\n";
      task->iterationIs(iteration);
    }

    // Jumps
    LocalNetwork::TaskList jumpProjectors = network->activeJumpProjectors();
    for (LocalNetwork::TaskList::iterator it = jumpProjectors.begin();
        it != jumpProjectors.end();
        ++it) {
      NamedTask::Ptr task = *it;
      log() << "Task: " << task->name() << "\n";
      task->iterationIs(iteration);
    }

    // Correction
    LocalNetwork::TaskList fullSliceTasks = network->activeFullTimeSlices();
    for (LocalNetwork::TaskList::iterator it = fullSliceTasks.begin();
        it != fullSliceTasks.end();
        ++it) {
      NamedTask::Ptr task = *it;
      log() << "Task: " << task->name() << "\n";
      task->iterationIs(iteration);
    }

    // Correction Synchronization 
    LocalNetwork::TaskList synchronizationTasks = network->activeCorrectionSyncs();
    for (LocalNetwork::TaskList::iterator it = synchronizationTasks.begin();
        it != synchronizationTasks.end();
        ++it) {
      NamedTask::Ptr task = *it;
      log() << "Task: " << task->name() << "\n";
      task->iterationIs(iteration);
    }

    // Update Seeds
    LocalNetwork::TaskList assemblyTasks = network->activeSeedAssemblers();
    for (LocalNetwork::TaskList::iterator it = assemblyTasks.begin();
        it != assemblyTasks.end();
        ++it) {
      NamedTask::Ptr task = *it;
      log() << "Task: " << task->name() << "\n";
      task->iterationIs(iteration);
    }

    // Local Fine Propagation
    LocalNetwork::TaskList halfSlices = network->activeHalfTimeSlices();
    for (LocalNetwork::TaskList::iterator it = halfSlices.begin();
        it != halfSlices.end();
        ++it) {
      NamedTask::Ptr task = *it;
      log() << "Task: " << task->name() << "\n";
      task->iterationIs(iteration);
    }

  }
}

void
ReducedLinearDriverImpl::solveCoarse(Communicator * timeComm) {
  RemoteCoarseServer::Ptr coarseServer;

  if (noForce_) {
    SeedInitializer::Ptr seedInitializer = buildSeedInitializer(true, NULL);
    coarseServer = RemoteSeedInitializerServer::New(timeComm, seedInitializer.ptr(), mapping_.ptr());

    log() << "\n";
    log() << "Remote coarse initialization\n";
  } else {
    DynamPropagator::Ptr coarsePropagator = buildCoarsePropagator(true, NULL);

    RemoteDynamPropagatorServer::Ptr propagatorServer = new RemoteDynamPropagatorServer(coarsePropagator.ptr(), timeComm);
    coarseServer = RemoteCoarseCorrectionServer::New(propagatorServer.ptr(), mapping_.ptr());

    mapping_->convergedSlicesInc(); // Since the coarse correction occurs at the second iteration (iteration 0 is after after iteration -1)
    
    log() << "\n";
    log() << "Remote coarse propagation\n";
  }

  coarseServer->statusIs(RemoteCoarseServer::BUSY);
}

DynamState
ReducedLinearDriverImpl::initialSeed() const {
  DynamState result = DynamState(vectorSize_);

  Vector & init_d = result.displacement();
  Vector & init_v = result.velocity();
  Vector init_a(vectorSize_);
  Vector init_vp(vectorSize_);
  domain()->initDispVeloc(init_d, init_v, init_a, init_vp);

  return result;
}

HalfTimeSlice::Manager::Ptr
ReducedLinearDriverImpl::buildHalfTimeSliceManager(GeneralizedAlphaParameter fineIntegrationParam,
                                                   PostProcessing::Manager * postProcessingMgr,
                                                   BasisCollector * collector) const {
  FineIntegratorManager::Ptr fineIntegratorMgr = LinearFineIntegratorManager<AffineGenAlphaIntegrator>::New(dynamOpsMgr_.ptr(), fineIntegrationParam);

  PropagatorManager::Ptr propagatorMgr =
    new PropagatorManager(
        collector,
        fineIntegratorMgr.ptr(), 
        postProcessingMgr,
        halfSliceRatio_,
        initialTime_);
  
  return HalfTimeSliceImpl::Manager::New(propagatorMgr.ptr());
}

PostProcessing::Manager::Ptr
ReducedLinearDriverImpl::buildPostProcessor(CpuRank localCpu) const {
  CpuRank myCpu();
  std::vector<int> localFileId;
  for (SliceMapping::SliceIterator it = mapping_->hostedSlice(localCpu); it; ++it) {
    localFileId.push_back((*it).value());
  }
  IncrementalPostProcessor::Ptr pitaPostProcessor = IncrementalPostProcessor::New(geoSource, localFileId.size(), &localFileId[0], probDesc()->getPostProcessor());
  typedef PostProcessing::IntegratorReactorImpl<IncrementalPostProcessor> LinearIntegratorReactor;
  return PostProcessing::Manager::New(LinearIntegratorReactor::Builder::New(pitaPostProcessor.ptr()).ptr());
}

BasisCollectorImpl::Ptr
ReducedLinearDriverImpl::buildBasisCollector() const {
  if (noForce_) {
    return BasisCollectorImpl::New();
  }

  return NonHomogeneousBasisCollectorImpl::New();
}

CorrectionTimeSlice::Manager::Ptr
ReducedLinearDriverImpl::buildCoarseCorrection(Communicator * coarseComm) const {
  if (noForce_) {
    return NULL;
  }

  DynamPropagator::Ptr coarsePropagator = buildCoarsePropagator(!remoteCoarse_, coarseComm);
  return LocalCorrectionTimeSlice::Manager::New(coarsePropagator.ptr());
}

LinearGenAlphaIntegrator::Ptr
ReducedLinearDriverImpl::buildCoarseIntegrator() const {
  return new HomogeneousGenAlphaIntegrator(dynamOpsMgr_.ptr(), GeneralizedAlphaParameter(coarseTimeStep_, coarseRhoInfinity_));
}

DynamPropagator::Ptr
ReducedLinearDriverImpl::buildCoarsePropagator(bool local, Communicator * coarseComm) const {
  if (local) {
    LinearGenAlphaIntegrator::Ptr coarseIntegrator = buildCoarseIntegrator();
    IntegratorPropagator::Ptr integratorPropagator = IntegratorPropagator::New(coarseIntegrator.ptr());
    integratorPropagator->timeStepCountIs(TimeStepCount(1));
    return integratorPropagator;
  }

  return RemoteDynamPropagatorProxy::New(vectorSize_, coarseComm, CpuRank(0));
}


SeedInitializer::Ptr
ReducedLinearDriverImpl::buildSeedInitializer(bool local, Communicator * timeComm) const {
  if (local) {
    LinearGenAlphaIntegrator::Ptr coarseIntegrator = buildCoarseIntegrator();
    coarseIntegrator->initialConditionIs(initialSeed(), initialTime_);
    return IntegratorSeedInitializer::New(coarseIntegrator.ptr(), TimeStepCount(1));
  }

  return RemoteSeedInitializerProxy::New(timeComm, vectorSize_);
}

} /* end namespace Hts */ } /* end namespace Pita */

// Global state
extern Communicator * structCom;
extern Domain * domain;

/* Entrypoint */
Pita::LinearDriver::Ptr
linearPitaDriverNew(SingleDomainDynamic * pbDesc) {
  return Pita::Hts::ReducedLinearDriverImpl::New(pbDesc, domain, &domain->solInfo(), structCom);
}

