#include "Fwk.h"
#include "Types.h"

#include "ReducedLinearDriverImpl.h"

#include "../DynamState.h"
#include <Problems.d/DynamDescr.h>

#include "../LinearGenAlphaIntegrator.h"
#include "../HomogeneousGenAlphaIntegrator.h"
#include "../IntegratorPropagator.h"
#include "../IntegratorSeedInitializer.h"

#include "SliceMapping.h"

#include "CorrectionNetworkImpl.h"

#include "../Seed.h"

#include "../RemoteStateMpiImpl.h"

#include "../IncrementalPostProcessor.h"
#include "LinearFineIntegratorManager.h"
#include "PropagatorManager.h"
#include "HalfTimeSliceImpl.h"

#include "ReducedFullTimeSliceImpl.h"
#include "../JumpProjector.h"
#include "../UpdatedSeedAssembler.h"
#include "LocalCorrectionTimeSlice.h"

#include <Comm.d/Communicator.h>
extern Communicator * structCom;


#include <Driver.d/Domain.h>
extern Domain * domain;

#include "LocalNetwork.h"

namespace Pita { namespace Hts {

ReducedLinearDriverImpl::ReducedLinearDriverImpl(SingleDomainDynamic<double> * pbDesc) :
  probDesc_(pbDesc)
{
  probDesc_->preProcess();
  
  /* Time domain */
  fineTimeStep_ = Seconds(domain->solInfo().dt);
  halfSliceRatio_ = TimeStepCount(domain->solInfo().Jratio / 2);
  sliceRatio_ = TimeStepCount(halfSliceRatio_.value() * 2);
  coarseTimeStep_ = fineTimeStep_ * sliceRatio_.value(); 
  initialTime_ = Seconds(domain->solInfo().initialTime);
  finalTime_ = Seconds(domain->solInfo().tmax);
  
  /* Time domain decomposition */
  numSlices_ = HalfSliceCount(static_cast<int>(ceil((finalTime_.value() - initialTime_.value()) / (halfSliceRatio_.value() * fineTimeStep_.value()))));
  fullTimeSlices_ = FullSliceCount((numSlices_.value() / 2) + (numSlices_.value() % 2));
  numSlices_ = HalfSliceCount(fullTimeSlices_.value() * 2);
  finalTime_ = fineTimeStep_ * Seconds(numSlices_.value() * halfSliceRatio_.value());

  /* Other parameters & options */ 
  maxActive_ = HalfSliceCount(domain->solInfo().numTSperCycleperCPU);
  lastIteration_ = IterationRank(domain->solInfo().kiter);
  projectorTolerance_ = domain->solInfo().pitaProjTol;
  noForce_ = domain->solInfo().NoForcePita;
 
  /* Available cpus and communication */ 
  remoteCoarse_ = false; //domain->solInfo().remoteCoarse && (timeCom_->numCPUs() > 1);
  timeCom_ = structCom;
  numCpus_ = CpuCount(timeCom_->numCPUs() - (remoteCoarse_ ? 1 : 0));
  myCpu_ = CpuRank(timeCom_->myID());
 
  /* Initial state */
  vectorSize_ = probDesc_->solVecInfo();
  initialSeed_ = DynamState(vectorSize_);
  Vector & init_d = initialSeed_.displacement();
  Vector & init_v = initialSeed_.velocity();
  Vector init_a(vectorSize_);
  Vector init_vp(vectorSize_);
  domain->initDispVeloc(init_d, init_v, init_a, init_vp);
}

void
ReducedLinearDriverImpl::solve() {
  double tic = getTime();
  
  /* Summarize problem and parameters */ 
  log() << "vectorSize = " << vectorSize_ << "\n";
  log() << "Slices = " << numSlices_ << ", MaxActive = " << maxActive_ << ", Cpus = " << numCpus_ << "\n";
  log() << "Num iter = " << lastIteration_ << "\n"; 
  log() << "dt = " << fineTimeStep_ << ", J/2 = " << halfSliceRatio_ << ", Dt = J*dt = " << coarseTimeStep_ << ", tf = " << finalTime_ << "\n";
 
  try {
    solveParallel();
  } catch(Fwk::Exception & e) {
    log() << e.what() << "\n";
  }
  
  double toc = getTime();

  log() << "\nTotal Solve Time = " << (toc - tic) / 1000.0 << " s\n";
  log() << "\nEnd Reduced HalfSlice Linear Pita\n";
}

void
ReducedLinearDriverImpl::solveParallel() {

  /* Time Integration */
  LinearDynamOps::Manager::Ptr dopsManager = LinearDynamOps::Manager::New(probDesc());
  double fineRhoInfinity = 1.0;

  GeneralizedAlphaParameter fineIntegrationParam(fineTimeStep_, fineRhoInfinity);
  LinearDynamOps::Ptr dynamOps = dopsManager->dynOpsNew(fineIntegrationParam);

  /* Postprocessing */ 
  std::vector<int> localFileId(numSlices_.value(), 0);
  for (int i = 0; i < numSlices_.value(); ++i) {
    localFileId[i] = i;
  }
  IncrementalPostProcessor::Ptr pitaPostProcessor = IncrementalPostProcessor::New(geoSource, localFileId.size(), &localFileId[0], probDesc_->getPostProcessor());
  typedef PostProcessing::IntegratorReactorImpl<IncrementalPostProcessor> LinearIntegratorReactor;
  PostProcessing::Manager::Ptr postProcessingMgr = PostProcessing::Manager::New(LinearIntegratorReactor::Builder::New(pitaPostProcessor.ptr()).ptr());

  /* Correction */
  SliceMapping::Ptr mapping = SliceMapping::New(fullTimeSlices_, numCpus_, maxActive_.value()); // TODO

  CorrectionNetworkImpl::Ptr correctionMgr =
    CorrectionNetworkImpl::New(
        vectorSize_,
        timeCom_,
        myCpu_,
        mapping.ptr(),
        dynamOps.ptr(),
        (noForce_ ? CorrectionNetworkImpl::HOMOGENEOUS : CorrectionNetworkImpl::NON_HOMOGENEOUS),
        projectorTolerance_); 

  /* Fine integrators */
  FineIntegratorManager::Ptr fineIntegratorMgr = LinearFineIntegratorManager<AffineGenAlphaIntegrator>::New(dopsManager.ptr(), fineIntegrationParam);

  PropagatorManager::Ptr propagatorMgr =
    new PropagatorManager(
        correctionMgr->collector(),
        fineIntegratorMgr.ptr(), 
        postProcessingMgr.ptr(),
        halfSliceRatio_,
        initialTime_);

  /* Tasks */
  HalfTimeSliceImpl::Manager::Ptr hsMgr = HalfTimeSliceImpl::Manager::New(propagatorMgr.ptr());
  JumpProjector::Manager::Ptr jpMgr = correctionMgr->jumpProjectorMgr();
  ReducedFullTimeSlice::Manager::Ptr fsMgr = correctionMgr->fullTimeSliceMgr();
  UpdatedSeedAssembler::Manager::Ptr usaMgr = correctionMgr->updatedSeedAssemblerMgr();

  /* Mpi communication manager */
  RemoteState::MpiManager::Ptr commMgr = RemoteState::MpiManager::New(timeCom_);
  commMgr->vectorSizeIs(vectorSize_);
 
  /* Coarse time-integration */ 
  double coarseRhoInfinity = 1.0;
  LinearGenAlphaIntegrator::Ptr coarseIntegrator = new HomogeneousGenAlphaIntegrator(dopsManager.ptr(), GeneralizedAlphaParameter(coarseTimeStep_, coarseRhoInfinity));
  coarseIntegrator->initialConditionIs(initialSeed_, initialTime_);
  DynamPropagator::Ptr coarsePropagator = IntegratorPropagator::New(coarseIntegrator.ptr());
 
  CorrectionTimeSlice::Manager::Ptr ctsMgr;
  if (!noForce_) {
    ctsMgr = LocalCorrectionTimeSlice::Manager::New(coarsePropagator.ptr());
  }
  
  /* Local tasks */
  LocalNetwork::Ptr network = new LocalNetwork(
      fullTimeSlices_,
      numCpus_,
      maxActive_,
      myCpu_,
      hsMgr.ptr(),
      jpMgr.ptr(),
      fsMgr.ptr(),
      usaMgr.ptr(),
      ctsMgr.ptr(),
      commMgr.ptr());

  log() << "\n";

  /* Initial seeds */
  IterationRank iteration;
 
  if (noForce_) { 
    SeedInitializer::Ptr seedInitializer = IntegratorSeedInitializer::New(coarseIntegrator.ptr(), TimeStepCount(1)); // TODO provided initial seeds
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
      seed->stateIs(seed->name() == "M0" ? initialSeed_ : DynamState(vectorSize_, 0.0));
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
      mapping->convergedSlicesInc(); // TODO remove one

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
    mapping->convergedSlicesInc(); // TODO remove one

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
solveCoarse() {
  // TODO
}

} /* end namespace Hts */ } /* end namespace Pita */

/* Entrypoint */
Pita::LinearDriver::Ptr
linearPitaDriverNew(SingleDomainDynamic<double> * pbDesc) {
  return Pita::Hts::ReducedLinearDriverImpl::New(pbDesc);
}
