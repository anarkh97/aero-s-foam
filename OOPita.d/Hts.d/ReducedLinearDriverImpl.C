#include "Fwk.h"
#include "Types.h"

#include "ReducedLinearDriverImpl.h"

#include "../DynamState.h"
#include <Problems.d/DynamDescr.h>

#include "../LinearGenAlphaIntegrator.h"
#include "../HomogeneousGenAlphaIntegrator.h"
#include "../IntegratorPropagator.h"
#include "../IntegratorSeedInitializer.h"

#include "StaticSliceStrategy.h"
#include "SliceMapping.h"

#include "LinearHalfSliceSchedule.h"
#include "CorrectionNetworkImpl.h"

#include "../Seed.h"

#include "../RemoteStateMpiImpl.h"

#include "../AffinePostProcessor.h"
#include "LinearFineIntegratorManager.h"
#include "PropagatorManager.h"
#include "HalfTimeSliceImpl.h"

#include "ReducedFullTimeSliceImpl.h"
#include "../JumpProjector.h"
#include "../UpdatedSeedAssembler.h"

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
  
  try {

  /* Summarize problem and parameters */ 
  log() << "vectorSize = " << vectorSize_ << "\n";
  log() << "Slices = " << numSlices_ << ", MaxActive = " << maxActive_ << ", Cpus = " << numCpus_ << "\n";
  log() << "Num iter = " << lastIteration_ << "\n"; 
  log() << "dt = " << fineTimeStep_ << ", J/2 = " << halfSliceRatio_ << ", Dt = J*dt = " << coarseTimeStep_ << ", tf = " << finalTime_ << "\n";

  /* Integration parameters */
  LinearDynamOps::Manager::Ptr dopsManager = LinearDynamOps::Manager::New(probDesc());
  double coarseRhoInfinity = 1.0;
  double fineRhoInfinity = 1.0;

  LinearGenAlphaIntegrator::Ptr coarseIntegrator = new HomogeneousGenAlphaIntegrator(dopsManager.ptr(), GeneralizedAlphaParameter(coarseTimeStep_, coarseRhoInfinity));
  coarseIntegrator->initialConditionIs(initialSeed_, initialTime_);
  SeedInitializer::Ptr seedInitializer = IntegratorSeedInitializer::New(coarseIntegrator.ptr(), TimeStepCount(1));
  DynamPropagator::Ptr coarsePropagator = IntegratorPropagator::New(coarseIntegrator.ptr());
  
  GeneralizedAlphaParameter fineIntegrationParam(fineTimeStep_, fineRhoInfinity);
  LinearDynamOps::Ptr dynamOps = dopsManager->dynOpsNew(fineIntegrationParam);
 
  /* Postprocessing */ 
  std::vector<int> localFileId(numSlices_.value(), 0);
  for (int i = 0; i < numSlices_.value(); ++i) {
    localFileId[i] = i;
  }
  AffinePostProcessor::Ptr pitaPostProcessor = AffinePostProcessor::New(geoSource, localFileId.size(), &localFileId[0], probDesc_->getPostProcessor());
  typedef PostProcessing::IntegratorReactorImpl<AffinePostProcessor> LinearIntegratorReactor;
  PostProcessing::Manager::Ptr postProcessingMgr = PostProcessing::Manager::New(LinearIntegratorReactor::Builder::New(pitaPostProcessor.ptr()).ptr());

  /* Correction */
  HalfSliceSchedule::Ptr schedule = LinearHalfSliceSchedule::New(numSlices_);
  SliceMapping::Ptr mapping = SliceMapping::New(fullTimeSlices_, numCpus_, maxActive_.value(), StaticSliceStrategy::New().ptr()); // TODO remove
  
  CorrectionNetworkImpl::Ptr correctionMgr =
    CorrectionNetworkImpl::New(
      vectorSize_,
      timeCom_,
      myCpu_,
      schedule.ptr(),
      mapping.ptr(),
      dynamOps.ptr(),
      CorrectionNetworkImpl::HOMOGENEOUS,
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

  /* Slice Managers */
  HalfTimeSliceImpl::Manager::Ptr hsMgr = HalfTimeSliceImpl::Manager::New(propagatorMgr.ptr());
  JumpProjector::Manager::Ptr jpMgr = correctionMgr->jumpProjectorMgr();
  ReducedFullTimeSlice::Manager::Ptr fsMgr = correctionMgr->fullTimeSliceMgr();
  UpdatedSeedAssembler::Manager::Ptr usaMgr = correctionMgr->updatedSeedAssemblerMgr();

  /* Mpi communication manager */
  RemoteState::MpiManager::Ptr commMgr = RemoteState::MpiManager::New(timeCom_);
  commMgr->vectorSizeIs(vectorSize_);

  LocalNetwork::Ptr network = new LocalNetwork(
      fullTimeSlices_,
      numCpus_,
      maxActive_,
      myCpu_,
      hsMgr.ptr(),
      jpMgr.ptr(),
      fsMgr.ptr(),
      usaMgr.ptr(),
      commMgr.ptr()); 

  log() << "\n";

  IterationRank iteration(0);
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

  log() << "\nIteration " << iteration << "\n";

  // HTS
  LocalNetwork::TaskList halfSlices = network->activeHalfTimeSlices();
  for (LocalNetwork::TaskList::iterator it = halfSlices.begin();
      it != halfSlices.end();
      ++it) {
    NamedTask::Ptr task = *it;
    log() << "Task: " << task->name() << "\n";
    task->iterationIs(iteration);
  }
 
  // Iteration loop
  for (iteration = IterationRank(1); iteration <= lastIteration_; iteration = iteration.next()) {

    log() << "\nIteration " << iteration << "\n";
    
    // Basis
    correctionMgr->buildProjection();
    size_t reducedBasisSize = correctionMgr->reducedBasisSize();
    log() << "ReducedBasisSize = " << reducedBasisSize << "\n";
    commMgr->reducedStateSizeIs(reducedBasisSize);

    // Next iteration 
    network->convergedSlicesInc();
    mapping->convergedSlicesInc(); // TODO remove

    // Propagated Seed Synchronization 
    LocalNetwork::TaskList leftSeedSyncs = network->activeLeftSeedSyncs();
    for (LocalNetwork::TaskList::iterator it = leftSeedSyncs.begin();
        it != leftSeedSyncs.end();
        ++it) {
      NamedTask::Ptr task = *it;
      log() << "Task: " << task->name() << "\n";
      task->iterationIs(iteration);
    }

    // Jump
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

    // Update
    LocalNetwork::TaskList assemblyTasks = network->activeSeedAssemblers();
    for (LocalNetwork::TaskList::iterator it = assemblyTasks.begin();
        it != assemblyTasks.end();
        ++it) {
      NamedTask::Ptr task = *it;
      log() << "Task: " << task->name() << "\n";
      task->iterationIs(iteration);
    }
  
    // HTS
    LocalNetwork::TaskList halfSlices = network->activeHalfTimeSlices();
    for (LocalNetwork::TaskList::iterator it = halfSlices.begin();
        it != halfSlices.end();
        ++it) {
      NamedTask::Ptr task = *it;
      log() << "Task: " << task->name() << "\n";
      task->iterationIs(iteration);
    }

  }

  /*LocalNetwork::HTSList activeList; //= network->activeHalfTimeSlices();

  log() << "Starting at slice " << network->firstActiveSlice() << "\n";
  activeList = network->activeHalfTimeSlices();
  for (LocalNetwork::HTSList::const_iterator it = activeList.begin(); it != activeList.end(); ++it) {
    log() << (*it)->rank() << (*it)->direction() << "\n";
  }

  log() << "Next iteration\n";
  network->convergedSlicesInc();

  log() << "Starting at slice " << network->firstActiveSlice() << "\n";
  activeList = network->activeHalfTimeSlices();
  for (LocalNetwork::HTSList::const_iterator it = activeList.begin(); it != activeList.end(); ++it) {
    log() << (*it)->rank() << (*it)->direction() << "\n";
  }*/


  /* Instantiate slices */
  /*for (HalfSliceRank rank(0); rank <= HalfSliceRank(0) + numSlices_; rank = rank + HalfSliceCount(1)) {
    Seed::Ptr mainSeed = seedMgr->instanceNew(String("MS_") + toString(rank)); 
    Seed::Ptr leftSeed = seedMgr->instanceNew(String("LPS_") + toString(rank)); 
    Seed::Ptr rightSeed = seedMgr->instanceNew(String("RPS_") + toString(rank));
    Seed::Ptr seedJump = seedMgr->instanceNew(String("SJ_") + toString(rank));

    ReducedSeed::Ptr jumpComp = reducedSeedMgr->instanceNew(String("JC_") + toString(rank)); 
    ReducedSeed::Ptr correctionComp = reducedSeedMgr->instanceNew(String("CC_") + toString(rank)); 
  }

  for (HalfSliceRank rank = HalfSliceRank(0); rank < HalfSliceRank(0) + numSlices_; rank = rank + HalfSliceCount(1)) {
    HalfTimeSliceImpl::Ptr forwardSlice = hsMgr->instanceNew(HalfSliceId(rank, HalfTimeSlice::FORWARD));
    forwardSlice->seedIs(seedMgr->instance(String("MS_") + toString(rank)));
    forwardSlice->propagatedSeedIs(seedMgr->instance(String("LPS_") + toString(rank + HalfSliceCount(1))));

    HalfTimeSliceImpl::Ptr backwardSlice = hsMgr->instanceNew(HalfSliceId(rank, HalfTimeSlice::BACKWARD));
    backwardSlice->seedIs(seedMgr->instance(String("MS_") + toString(rank + HalfSliceCount(1))));
    backwardSlice->propagatedSeedIs(seedMgr->instance(String("RPS_") + toString(rank)));
    
    JumpProjector::Ptr jumpProjector = jpMgr->instanceNew(toString(rank));
    jumpProjector->predictedSeedIs(seedMgr->instance(String("RPS_") + toString(rank)));
    jumpProjector->actualSeedIs(seedMgr->instance(String("LPS_") + toString(rank)));
    jumpProjector->seedJumpIs(seedMgr->instance(String("SJ_") + toString(rank)));
    jumpProjector->reducedSeedJumpIs(reducedSeedMgr->instance(String("JC_") + toString(rank)));

    UpdatedSeedAssembler::Ptr seedAssembler = usMgr->instanceNew(toString(rank));
    seedAssembler->updatedSeedIs(seedMgr->instance(String("MS_") + toString(rank)));
    seedAssembler->propagatedSeedIs(seedMgr->instance(String("LPS_") + toString(rank)));
    seedAssembler->correctionComponentsIs(reducedSeedMgr->instance(String("CC_") + toString(rank)));
  }

  for (HalfSliceRank rank = HalfSliceRank(0); rank < HalfSliceRank(0) + numSlices_ - HalfSliceCount(1); rank = rank + HalfSliceCount(1)) {
    ReducedFullTimeSlice::Ptr fullSlice = fsMgr->instanceNew(rank);
    fullSlice->jumpIs(reducedSeedMgr->instance(String("JC_") + toString(rank)));
    fullSlice->correctionIs(reducedSeedMgr->instance(String("CC_") + toString(rank)));
    fullSlice->nextCorrectionIs(reducedSeedMgr->instance(String("CC_") + toString(rank + HalfSliceCount(2))));
  }*/

  /* Classify slices */
  /*struct IterationData {
    std::deque<HalfTimeSliceImpl::Ptr> halfSlices;
    std::deque<JumpProjector::Ptr> jumpProjectors;
    std::deque<ReducedFullTimeSlice::Ptr> fullSlices;
    std::deque<UpdatedSeedAssembler::Ptr> seedAssemblers;
  };
  IterationData data[2];

  std::deque<Seed::Ptr> primalMainSeeds;
  primalMainSeeds.push_back(seedMgr->instance("MS_0"));

  for (FullSliceRank rank(0); rank < FullSliceRank(0) + fullTimeSlices_; rank = rank + FullSliceCount(1)) {
    HalfSliceRank headRank = HalfSliceRank(rank.value() * 2);
    HalfSliceRank tailRank = headRank + HalfSliceCount(1);

    data[0].halfSlices.push_back(hsMgr->instance(HalfSliceId(headRank, HalfTimeSlice::FORWARD)));
    data[0].halfSlices.push_back(hsMgr->instance(HalfSliceId(tailRank, HalfTimeSlice::BACKWARD)));

    data[0].jumpProjectors.push_back(jpMgr->instance(toString(tailRank)));
    ReducedFullTimeSlice::Ptr fs = fsMgr->instance(tailRank);
    if (fs) data[0].fullSlices.push_back(fs);
    data[0].seedAssemblers.push_back(usMgr->instance(toString(tailRank)));

    data[1].halfSlices.push_back(hsMgr->instance(HalfSliceId(headRank, HalfTimeSlice::BACKWARD)));
    data[1].halfSlices.push_back(hsMgr->instance(HalfSliceId(tailRank, HalfTimeSlice::FORWARD)));
    
    data[1].jumpProjectors.push_back(jpMgr->instance(toString(headRank)));
    fs = fsMgr->instance(headRank);
    if (fs) data[1].fullSlices.push_back(fs);
    data[1].seedAssemblers.push_back(usMgr->instance(toString(headRank)));

    primalMainSeeds.push_back(seedMgr->instance(String("MS_" + toString(tailRank + HalfSliceCount(1)))));
  }*/

  /* Do the stuff */
  // Initialize
  /*for (std::deque<Seed::Ptr>::iterator it = primalMainSeeds.begin();
       it != primalMainSeeds.end();
       ++it) {
    FullSliceRank rank(std::distance(primalMainSeeds.begin(), it));
    //log() << "Initial Seed " << rank << " " << (*it)->name() << "\n";
    DynamState initSeed = seedInitializer->initialSeed(rank);
    (*it)->stateIs(initSeed);
    (*it)->iterationIs(IterationRank(0));
    (*it)->statusIs(Seed::ACTIVE);
  }

  // Loop
  //for (IterationRank iteration(0); iteration <= lastIteration_; iteration = iteration.next()) {
  for (IterationRank iteration(0); iteration < lastIteration_; iteration = iteration.next()) {
    int dataIndex = iteration.value() % 2;
    log() << "\nIteration # " << iteration << "\n";
    IterationData & iterData = data[dataIndex];

    for (std::deque<HalfTimeSliceImpl::Ptr>::iterator it = iterData.halfSlices.begin();
        it != iterData.halfSlices.end();
        ++it) {
      //log() << "HTS ";
      //log() << (*it)->rank() << ((*it)->direction() == HalfTimeSlice::FORWARD ? 'F' : 'B') << "\n";
      (*it)->propagateSeed();
    }

    // Basis
    correctionMgr->buildProjection();

    // Correction
    for (std::deque<JumpProjector::Ptr>::iterator it = iterData.jumpProjectors.begin();
        it != iterData.jumpProjectors.end();
        ++it) {
      (*it)->iterationIs(iteration);
    }

    ReducedSeed::Ptr firstCorrection = reducedSeedMgr->instance(String("CC_") + toString(mapping->firstActiveSlice() + HalfSliceCount(1)));
    firstCorrection->stateIs(Vector(iterData.fullSlices.front()->jump()->state().size(), 0.0));
    log() << "First correction is " << firstCorrection->name() << "\n";

    for (std::deque<ReducedFullTimeSlice::Ptr>::iterator it = iterData.fullSlices.begin();
        it != iterData.fullSlices.end();
        ++it) {
      (*it)->iterationIs(iteration);
    }

    for (std::deque<UpdatedSeedAssembler::Ptr>::iterator it = iterData.seedAssemblers.begin();
        it != iterData.seedAssemblers.end();
        ++it) {
      (*it)->doAssembly();
    }

    // Next iteration
    mapping->convergedSlicesInc();

    data[0].halfSlices.pop_front();
    data[1].halfSlices.pop_front();

    data[1 - dataIndex].jumpProjectors.pop_front();
    data[1 - dataIndex].fullSlices.pop_front();
    data[1 - dataIndex].seedAssemblers.pop_front(); 
  }*/

  } catch(Fwk::Exception & e) {
    log() << e.what() << "\n";
  }

  /* End */ 
  double toc = getTime();

  log() << "Total Solve Time = " << (toc - tic) / 1000.0 << " s\n";
  log() << "\nEnd Reduced HalfSlice Linear Pita\n";
}

} /* end namespace Hts */ } /* end namespace Pita */

/* Entrypoint */
Pita::LinearDriver::Ptr
linearPitaDriverNew(SingleDomainDynamic<double> * pbDesc) {
  return Pita::Hts::ReducedLinearDriverImpl::New(pbDesc);
}
