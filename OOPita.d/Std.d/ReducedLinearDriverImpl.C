#include "Fwk.h"
#include "Types.h"

#include "ReducedLinearDriverImpl.h"

#include "../DynamState.h"
#include "../LinearDynamOps.h"
#include "../DynamStateOps.h"

#include <Problems.d/DynamDescr.h>
#include <Utils.d/SolverInfo.h>
#include <Driver.d/Domain.h>

#include "SliceMapping.h"

#include "../Seed.h"
#include "IncrementalPropagation.h"

#include "../PostProcessingManager.h"
#include "../IncrementalPostProcessor.h"

#include "../HomogeneousGenAlphaIntegrator.h"
#include "LinearPropagatorManager.h"

#include "LinearProjectionNetwork.h"
#include "../NearSymmetricSolver.h"

#include "../JumpBuilder.h"
#include "../JumpProjection.h"
#include "../ReducedCorrectionPropagatorImpl.h"
#include "../FullCorrectionPropagatorImpl.h"
#include "../UpdatedSeedAssemblerImpl.h"
#include "JumpConvergenceEvaluator.h"

#include "../IntegratorSeedInitializer.h"
#include "RemoteSeedInitializerServer.h"
#include "../RemoteSeedInitializerProxy.h"
#include "../UserProvidedSeedInitializer.h"

#include "HomogeneousTaskManager.h"
#include "NonHomogeneousTaskManager.h"
#include "../TimedExecution.h"

#include "../RemoteStateMpiImpl.h"

#include "../CommSplitter.h"

#include "../IntegratorPropagator.h"
#include "../RemoteDynamPropagatorProxy.h"
#include "RemoteCoarseCorrectionServer.h"

#include <Timers.d/GetTime.h>

namespace Pita { namespace Std {

ReducedLinearDriverImpl::ReducedLinearDriverImpl(SingleDomainDynamic * pbDesc,
                                                 GeoSource * geoSource,
                                                 Domain * domain,
                                                 SolverInfo * solverInfo,
                                                 Communicator * baseComm) :
  LinearDriverImpl(pbDesc, geoSource, domain, solverInfo, baseComm)
{}

// Main routine
void
ReducedLinearDriverImpl::solve() {
  log() << "Begin Linear Pita\n";
  double tic = getTime(); // Total time

  preprocess();

  /* Summarize problem and parameters */
  log() << "\n"; 
  log() << "Slices = " << mapping_->totalSlices() << ", MaxActive = " << mapping_->maxWorkload() << ", Cpus = " << mapping_->availableCpus() << "\n";
  log() << "Iteration count = " << lastIteration_ << "\n"; 
  log() << "dt = " << fineTimeStep_ << ", J = " << sliceRatio_ << ", Dt = J*dt = " << coarseTimeStep_ << ", Tf = Slices*Dt = " << finalTime_ << "\n";
  if (noForce_) { log() << "No external force\n"; }
  if (remoteCoarse_) { log() << "Remote coarse time-integration\n"; }
  if (userProvidedSeeds_) { log() << "Reading user-provided initial seed information\n"; }
  log() << "VectorSize = " << vectorSize_ << " dofs\n";
  log() << "Projector tol = " << projectorTolerance_ << "\n";

  // Determine process task
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
  log() << "\n" << "Total time = " << (toc - tic) / 1000.0 << " s\n";
  log() << "\n" << "End Linear Pita\n";
}

void
ReducedLinearDriverImpl::preprocess() {
  double tic = getTime();
  
  probDesc()->preProcess();
  
  // Space-domain 
  vectorSize_ = probDesc()->solVecInfo();
  dynamOpsMgr_= LinearDynamOps::Manager::New(probDesc());

  // Time-domain 
  fineTimeStep_ = Seconds(solverInfo()->getTimeStep());
  sliceRatio_ = TimeStepCount(solverInfo()->pitaTimeGridRatio);
  coarseTimeStep_ = fineTimeStep_ * sliceRatio_.value(); 
  initialTime_ = Seconds(solverInfo()->initialTime);
  finalTime_ = Seconds(solverInfo()->tmax);
 
  SliceCount numSlices(static_cast<int>(ceil((finalTime_.value() - initialTime_.value()) / (sliceRatio_.value() * fineTimeStep_.value()))));
  finalTime_ = fineTimeStep_ * Seconds(numSlices.value() * sliceRatio_.value()); // To have a whole number of full time-slices 
  
  // Main options 
  noForce_ = solverInfo()->pitaNoForce;
  userProvidedSeeds_ = solverInfo()->pitaReadInitSeed && noForce_;
  remoteCoarse_ = solverInfo()->pitaRemoteCoarse && (baseComm()->numCPUs() > 1) && !userProvidedSeeds_;

  // Load balancing 
  CpuCount numCpus(baseComm()->numCPUs() - (remoteCoarse_ ? 1 : 0));
  SliceCount maxActive(solverInfo()->pitaProcessWorkloadMax);
  mapping_ = SliceMapping::New(numSlices, numCpus, maxActive);

  // Other parameters 
  lastIteration_ = IterationRank(solverInfo()->pitaMainIterMax);
  projectorTolerance_ = solverInfo()->pitaProjTol;
  coarseRhoInfinity_ = 1.0; // TODO Could be set in input file
 
  double toc = getTime();
  log() << "\n";
  log() << "Total preprocessing time = " << (toc - tic) / 1000.0 << " s\n";
}


void
ReducedLinearDriverImpl::solveParallel(Communicator * timeComm, Communicator * coarseComm) {
  log() << "\n";
  log() << "Parallel solver initialization\n";

  double tic = getTime();
   
  // Local process identity
  CpuRank myCpu(timeComm->myID());

  // Linear operators
  const double fineRhoInfinity = 1.0; // No damping for the time being
  GeneralizedAlphaParameter fineIntegrationParam(fineTimeStep_, fineRhoInfinity);
  LinearDynamOps::Ptr dynOps = dynamOpsMgr_->dynOpsNew(fineIntegrationParam); // Most expensive step 

  // Projection network
  RankDeficientSolver::Ptr normalMatrixSolver = NearSymmetricSolver::New(projectorTolerance_);
  LinearProjectionNetwork::Ptr projectionMgr = new LinearProjectionNetwork(mapping_.ptr(),
                                                                           timeComm,
                                                                           dynOps.ptr(),
                                                                           vectorSize_,
                                                                           normalMatrixSolver.ptr());
  // Correction operations
  JumpProjection::Manager::Ptr jumpProjMgr = JumpProjection::Manager::New(projectionMgr->projectionBasis());
  CorrectionPropagator<Vector>::Manager::Ptr corrPropMgr = ReducedCorrectionPropagatorImpl::Manager::New(
      projectionMgr->reprojectionMatrix(), projectionMgr->normalMatrixSolver());
  UpdatedSeedAssembler::Manager::Ptr seedUpMgr = UpdatedSeedAssemblerImpl::Manager::New(projectionMgr->propagatedBasis());
  
  // Fine-grid time integration 
  AffineGenAlphaIntegrator::Ptr fineIntegrator = new AffineGenAlphaIntegrator(dynamOpsMgr_.ptr(), fineIntegrationParam);
  PostProcessing::Manager::Ptr postProcessingMgr = buildPostProcessor(myCpu);
  AffineBasisCollector::Ptr collector = projectionMgr->collector();
  LinearPropagatorManager::Ptr finePropagatorManager = LinearPropagatorManager::New(
      fineIntegrator.ptr(), postProcessingMgr.ptr(), collector.ptr(), sliceRatio_, initialTime_);

  // Point-to-point communication
  RemoteState::MpiManager::Ptr commMgr = RemoteState::MpiManager::New(timeComm, vectorSize_);

  // Jump-based convergence policy
  bool useCvg = true; // TODO Parser option
  JumpConvergenceEvaluator::Ptr jumpCvgEval;
  if (useCvg) {
    const int schemeOrder = 2;
    const double jumpCvgRatio = std::pow(static_cast<double>(sliceRatio_.value()), schemeOrder);
    jumpCvgEval = new AccumulatedJumpConvergenceEvaluator(jumpCvgRatio, dynOps.ptr(), mapping_.ptr(), timeComm); 
  } else {
    jumpCvgEval = new TrivialConvergenceEvaluator(mapping_.ptr());
  }

  // Tasks
  TaskManager::Ptr taskMgr;
  if (noForce_) {
    // Initial seed information on coarse time-grid
    SeedInitializer::Ptr seedInitializer = buildSeedInitializer(coarseComm);
    taskMgr = new HomogeneousTaskManager(mapping_.ptr(),
                                         myCpu,
                                         commMgr.ptr(),
                                         seedInitializer.ptr(),
                                         finePropagatorManager.ptr(),
                                         projectionMgr.ptr(),
                                         jumpProjMgr.ptr(),
                                         corrPropMgr.ptr(),
                                         seedUpMgr.ptr(),
                                         jumpCvgEval.ptr());
  } else {
    // Coarse time-grid propagator
    CorrectionPropagator<DynamState>::Manager::Ptr fullCorrPropMgr = buildCoarseCorrection(coarseComm);
    taskMgr = new NonHomogeneousTaskManager(mapping_.ptr(),
                                            myCpu,
                                            commMgr.ptr(),
                                            initialSeed(),
                                            finePropagatorManager.ptr(),
                                            projectionMgr.ptr(),
                                            jumpProjMgr.ptr(),
                                            corrPropMgr.ptr(),
                                            fullCorrPropMgr.ptr(),
                                            seedUpMgr.ptr(),
                                            jumpCvgEval.ptr());
  }

  TimedExecution::Ptr execution = TimedExecution::New(taskMgr.ptr());
  
  double toc = getTime();
  log() << "\n";
  log() << "Total initialization time = " << (toc - tic) / 1000.0 << " s\n";
  tic = toc;
  
  execution->targetIterationIs(lastIteration_);

  toc = getTime();
  log() << "\n";
  log() << "Total solve time = " << (toc - tic) / 1000.0 << " s\n";
}


void
ReducedLinearDriverImpl::solveCoarse(Communicator * clientComm) {
  double tic = getTime();
  RemoteCoarseServer::Ptr coarseServer;

  if (noForce_) {
    SeedInitializer::Ptr seedInitializer = buildSeedInitializer();
    coarseServer = RemoteSeedInitializerServer::New(clientComm, seedInitializer.ptr(), mapping_.ptr());

    log() << "\n";
    log() << "Remote coarse initialization\n";
  } else {
    DynamPropagator::Ptr coarsePropagator = buildCoarsePropagator();

    RemoteDynamPropagatorServer::Ptr propagatorServer = new RemoteDynamPropagatorServer(coarsePropagator.ptr(), clientComm);
    coarseServer = RemoteCoarseCorrectionServer::New(propagatorServer.ptr(), mapping_.ptr());

    mapping_->convergedSlicesInc(); // Since the coarse correction occurs at the second iteration (iteration 0 after iteration -1)
    
    log() << "\n";
    log() << "Remote coarse propagation\n";
  }
  double toc = getTime();
  log() << "\n";
  log() << "Total initialization time = " << (toc - tic) / 1000.0 << " s\n";
  tic = getTime();
  
  coarseServer->statusIs(RemoteCoarseServer::BUSY);
  
  toc = getTime();
  log() << "\n";
  log() << "Total solve time = " << (toc - tic) / 1000.0 << " s\n";
}


PostProcessing::Manager::Ptr
ReducedLinearDriverImpl::buildPostProcessor(CpuRank localCpu) const {
  std::vector<int> localFileId;
  for (SliceMapping::SliceIterator it = mapping_->hostedSlice(localCpu); it; ++it) {
    localFileId.push_back((*it).value());
  }
  IncrementalPostProcessor::Ptr pitaPostProcessor = IncrementalPostProcessor::New(geoSource(), localFileId.size(), &localFileId[0], probDesc()->getPostProcessor());
  typedef PostProcessing::IntegratorReactorImpl<IncrementalPostProcessor> LinearIntegratorReactor;
  return PostProcessing::Manager::New(LinearIntegratorReactor::Builder::New(pitaPostProcessor.ptr()).ptr());
}


SeedInitializer::Ptr
ReducedLinearDriverImpl::buildSeedInitializer(Communicator * clientComm) const {
  if (userProvidedSeeds_) {
    return UserProvidedSeedInitializer::New(vectorSize_, geoSource(), domain());
  }

  if (clientComm == NULL) {
    // Local time-integration
    LinearGenAlphaIntegrator::Ptr coarseIntegrator = buildCoarseIntegrator();
    coarseIntegrator->initialConditionIs(initialSeed(), initialTime_);
    return IntegratorSeedInitializer::New(coarseIntegrator.ptr(), TimeStepCount(1));
  }

  return RemoteSeedInitializerProxy::New(clientComm, vectorSize_);
}

LinearGenAlphaIntegrator::Ptr
ReducedLinearDriverImpl::buildCoarseIntegrator() const {
  return new HomogeneousGenAlphaIntegrator(dynamOpsMgr_.ptr(), GeneralizedAlphaParameter(coarseTimeStep_, coarseRhoInfinity_));
}


CorrectionPropagator<DynamState>::Manager::Ptr
ReducedLinearDriverImpl::buildCoarseCorrection(Communicator * coarseComm) const {
  DynamPropagator::Ptr coarsePropagator = buildCoarsePropagator(coarseComm);
  return FullCorrectionPropagatorImpl::Manager::New(coarsePropagator.ptr());
}

DynamPropagator::Ptr
ReducedLinearDriverImpl::buildCoarsePropagator(Communicator * coarseComm) const {
  if (!coarseComm) {
    LinearGenAlphaIntegrator::Ptr coarseIntegrator = buildCoarseIntegrator();
    IntegratorPropagator::Ptr localPropagator = IntegratorPropagator::New(coarseIntegrator.ptr());
    localPropagator->timeStepCountIs(TimeStepCount(1));
    return localPropagator;
  }

  return RemoteDynamPropagatorProxy::New(vectorSize_, coarseComm, CpuRank(0));
}

} /* end namespace Std */ } /* end namespace Pita */

// Global state
extern GeoSource * geoSource;
extern Communicator * structCom;
extern Domain * domain;

/* Entrypoint */
Pita::LinearDriver::Ptr
linearPitaDriverNew(SingleDomainDynamic * pbDesc) {
  return Pita::Std::ReducedLinearDriverImpl::New(pbDesc, geoSource, domain, &domain->solInfo(), structCom);
}
