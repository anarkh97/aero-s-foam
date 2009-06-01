#include "LinearDriverImpl.h"

#include "Activity.h"

#include "HalfSliceNetwork.h"
#include "LinearHalfSliceSchedule.h"

#include "SliceMapping.h"
#include "StaticSliceStrategy.h"

#include "HalfTimeSliceImpl.h"
#include "FullTimeSliceImpl.h"

#include "HalfSliceCorrectionNetworkImpl.h"
#include "NonHomogeneousBasisCollectorImpl.h"
#include "HalfSlicePropagatorManager.h"
#include "LinearFineIntegratorManager.h"

#include "PostProcessingManager.h"

#include "SimpleSeedInitializer.h"
#include "IntegratorSeedInitializer.h"
#include "RemoteSeedInitializer.h"
#include "RemoteSeedInitializerProxy.h"

#include <Problems.d/DynamDescr.h> 

#include <Comm.d/Communicator.h>
extern Communicator * structCom;

#include <Driver.d/Domain.h>
extern Domain * domain;

#include "LinearGenAlphaIntegrator.h"
#include "HomogeneousGenAlphaIntegrator.h"

#include "AffinePostProcessor.h"

#include "Seed.h"
#include "ScheduledRemoteSeedWriter.h"
#include "ScheduledRemoteSeedReader.h"

#include "RemoteCoarseCorrectionServer.h"
#include "RemoteDynamPropagator.h"

#include "Log.h"
#include <Timers.d/GetTime.h>

namespace Pita {

// TODO HACK
using namespace Hts;

LinearDriverImpl::LinearDriverImpl(SingleDomainDynamic<double> * pbDesc) :
  probDesc_(pbDesc)
{
  probDesc_->preProcess();
}

void
LinearDriverImpl::solve() {

  double tic = getTime();

  Activity::Manager::Ptr activityMgr = activityManagerInstance();
  
  size_t vectorSize = probDesc_->solVecInfo();
  
  Seconds fineTimeStep(domain->solInfo().dt);
  TimeStepCount halfSliceRatio(domain->solInfo().Jratio / 2);
  TimeStepCount sliceRatio(halfSliceRatio.value() * 2);
  Seconds coarseTimeStep = fineTimeStep * sliceRatio.value(); 

  Seconds initialTime(domain->solInfo().initialTime);
  Seconds finalTime(domain->solInfo().tmax);
  HalfSliceCount numSlices(static_cast<int>(ceil((finalTime.value() - initialTime.value()) / (halfSliceRatio.value() * fineTimeStep.value()))));
  FullSliceCount fullTimeSlices = FullSliceCount((numSlices.value() / 2) + (numSlices.value() % 2));
  
  numSlices = HalfSliceCount(fullTimeSlices.value() * 2);
  finalTime = Seconds(numSlices.value() * halfSliceRatio.value() * fineTimeStep.value());

  HalfSliceSchedule::Ptr schedule = LinearHalfSliceSchedule::New(numSlices);
  
  IterationRank lastIteration(domain->solInfo().kiter);
  
  Communicator * timeCom = structCom;
  
  HalfSliceCount maxActive(domain->solInfo().numTSperCycleperCPU);

  bool remoteCoarse = domain->solInfo().remoteCoarse && (timeCom->numCPUs() > 1);
  CpuCount numCpus(timeCom->numCPUs() - (remoteCoarse ? 1 : 0));
  CpuRank myCpu(timeCom->myID());

  log() << "vectorSize = " << vectorSize << "\n";
  log() << "Slices = " << numSlices << ", MaxActive = " << maxActive << ", Cpus = " << numCpus << "\n";
  log() << "Num iter = " << lastIteration << "\n"; 
  log() << "dt = " << fineTimeStep << ", J/2 = " << halfSliceRatio << ", Dt = J*dt = " << coarseTimeStep << ", tf = " << finalTime << "\n";
  
  StaticSliceStrategy::Ptr strategy = StaticSliceStrategy::New();
  SliceMapping::Ptr mapping = SliceMapping::New(fullTimeSlices, numCpus, maxActive.value(), strategy.ptr());

  // No damping
  double rho_infinity_coarse = 1.0;
  double rho_infinity_fine = 1.0;
  
  LinearDynamOps::Manager::Ptr dopsManager = LinearDynamOps::Manager::New(probDesc());
  
  DynamState initialSeed(vectorSize);
  Vector & init_d = initialSeed.displacement();
  Vector & init_v = initialSeed.velocity();
  Vector init_a(vectorSize);
  Vector init_vp(vectorSize);
  domain->initDispVeloc(init_d, init_v, init_a, init_vp);
  
  HalfSliceCorrectionNetworkImpl::Strategy correctionStrategy = domain->solInfo().NoForcePita ?
    HalfSliceCorrectionNetworkImpl::HOMOGENEOUS :
    HalfSliceCorrectionNetworkImpl::NON_HOMOGENEOUS;

  SeedInitializer::Ptr seedInitializer;
  DynamPropagator::Ptr coarsePropagator;

  if (remoteCoarse) {
    int originId = timeCom->myID();
    int originCpus = timeCom->numCPUs();

    MPI_Comm * originCom = timeCom->getCommunicator();

    MPI_Comm splitCom;
    int color = (originId != 0) ? 1 : 0;
    MPI_Comm_split(*originCom, color, originId, &splitCom);

    int splitCpus, splitId;
    MPI_Comm_size(splitCom, &splitCpus);
    MPI_Comm_rank(splitCom, &splitId); 

    MPI_Comm interCom;
    int remoteLeaderId = (color == 0) ? 1 : 0;
    MPI_Intercomm_create(splitCom, 0, *originCom, remoteLeaderId, 0, &interCom);

    if (color == 0) {
      LinearGenAlphaIntegrator::Ptr coarseIntegrator = new HomogeneousGenAlphaIntegrator(dopsManager.ptr(), GeneralizedAlphaParameter(coarseTimeStep , rho_infinity_coarse));
      
      double tac = getTime();
      log() << "Total Init Time = " << (tac - tic) / 1000.0 << " s\n";
     
      Communicator * seedInitCom = new Communicator(interCom, stderr); // TODO delete
      
      if (correctionStrategy == HalfSliceCorrectionNetworkImpl::HOMOGENEOUS) {
        coarseIntegrator->initialConditionIs(initialSeed, initialTime);
        SeedInitializer::Ptr seedInitializer = IntegratorSeedInitializer::New(coarseIntegrator.ptr(), TimeStepCount(1));

        RemoteSeedInitializer::Ptr remoteSeedInitializer = RemoteSeedInitializer::New(seedInitCom, seedInitializer.ptr(), mapping.ptr());

        remoteSeedInitializer->statusIs(RemoteSeedInitializer::BUSY);

      } else {
        IntegratorPropagator::Ptr integratorPropagator = IntegratorPropagator::New(coarseIntegrator.ptr());
        integratorPropagator->timeStepCountIs(TimeStepCount(1));
        coarsePropagator = integratorPropagator;
        RemoteDynamPropagatorServer::Ptr propagatorServer = new RemoteDynamPropagatorServer(coarsePropagator.ptr(), seedInitCom);
        RemoteCoarseCorrectionServer::Ptr correctionServer = new RemoteCoarseCorrectionServer(propagatorServer.ptr(), mapping.ptr(), schedule->correction());
        correctionServer->statusIs(RemoteCoarseCorrectionServer::ACTIVE);
        
        activityMgr->targetPhaseIs(lastIteration, PhaseRank(0));
      }
        
      double toc = getTime();
      log() << "Total Solve Time = " << (toc - tic) / 1000.0 << " s\n";

      log() << "End LinPita\n";
      return;

    } else {
      timeCom = new Communicator(splitCom, stderr); // TODO delete

      myCpu = CpuRank(timeCom->myID());
      Communicator * seedInitCom = new Communicator(interCom, stderr); // TODO delete

      if (correctionStrategy == HalfSliceCorrectionNetworkImpl::HOMOGENEOUS) {
        seedInitializer = RemoteSeedInitializerProxy::New(seedInitCom, vectorSize);
      } else {
        seedInitializer = SimpleSeedInitializer::New(initialSeed);
        coarsePropagator = RemoteDynamPropagator::New(vectorSize, seedInitCom, CpuRank(0));
      }
    }
  } else {
    LinearGenAlphaIntegrator::Ptr coarseIntegrator = new HomogeneousGenAlphaIntegrator(dopsManager.ptr(), GeneralizedAlphaParameter(coarseTimeStep , rho_infinity_coarse));
    coarseIntegrator->initialConditionIs(initialSeed, initialTime);
    if (correctionStrategy == HalfSliceCorrectionNetworkImpl::HOMOGENEOUS) {
      seedInitializer = IntegratorSeedInitializer::New(coarseIntegrator.ptr(), TimeStepCount(1)); // TODO other coarse initialization
    } else {
      seedInitializer = SimpleSeedInitializer::New(initialSeed);
    }
    coarsePropagator = IntegratorPropagator::New(coarseIntegrator.ptr());
  } 

  RemoteHalfSliceSurrogate::Ptr surrogate = RemoteHalfSliceSurrogate::New(mapping.ptr());

  Seed::Manager::Ptr seedMgr = Seed::Manager::New();
  ScheduledRemoteSeedWriter::Manager::Ptr writerMgr = ScheduledRemoteSeedWriter::Manager::New(timeCom, surrogate.ptr());
  ScheduledRemoteSeedReader::Manager<Hs::CommId>::Ptr readerMgr = ScheduledRemoteSeedReader::Manager<Hs::CommId>::New(timeCom);

  GeneralizedAlphaParameter integrationParam(fineTimeStep, rho_infinity_fine);
  LinearDynamOps::Ptr dynamOps = dopsManager->dynOpsNew(integrationParam);

  HalfSliceCorrectionNetworkImpl::Ptr correctionMgr =
    HalfSliceCorrectionNetworkImpl::New(
      vectorSize,
      timeCom,
      myCpu,
      schedule.ptr(),
      mapping.ptr(),
      dynamOps.ptr(),
      correctionStrategy); 

  std::vector<int> localFileId;
  for (SliceMapping::SliceIdIterator it = mapping->hostedSlice(myCpu, HalfSliceRank(0), HalfSliceRank(0) + mapping->totalSlices()); it; ++it) {
    if (it->type() == Hs::FORWARD_HALF_SLICE || it->type() == Hs::BACKWARD_HALF_SLICE) {
      localFileId.push_back(it->rank().value());
    }
  }
  std::sort(localFileId.begin(), localFileId.end());
  localFileId.erase(std::unique(localFileId.begin(), localFileId.end()), localFileId.end());

  // HACK 
  LinearFineIntegratorManager<AffineGenAlphaIntegrator>::Ptr fineIntegratorMgr = LinearFineIntegratorManager<AffineGenAlphaIntegrator>::New(dopsManager.ptr(), integrationParam);
  AffinePostProcessor::Ptr pitaPostProcessor = AffinePostProcessor::New(geoSource, localFileId.size(), &localFileId[0], probDesc_->getPostProcessor());
  typedef PostProcessing::IntegratorReactorImpl<AffinePostProcessor> LinearIntegratorReactor;
  PostProcessing::Manager::Ptr postProcessingMgr = PostProcessing::Manager::New(LinearIntegratorReactor::Builder::New(pitaPostProcessor.ptr()).ptr());

  HalfSlicePropagatorManager::Ptr propagatorMgr =
    new HalfSlicePropagatorManager(
      correctionMgr->collector(),
      fineIntegratorMgr.ptr(), 
      postProcessingMgr.ptr(),
      halfSliceRatio,
      initialTime);

  HalfTimeSliceImpl::Manager::Ptr hsMgr = HalfTimeSliceImpl::Manager::New(propagatorMgr.ptr());
  FullTimeSliceHeadImpl::Manager::Ptr fshMgr = FullTimeSliceHeadImpl::Manager::New(correctionMgr->reductorMgr(), coarsePropagator.ptr(), seedMgr.ptr(), timeCom, myCpu);
  FullTimeSliceTailImpl::Manager::Ptr fstMgr = FullTimeSliceTailImpl::Manager::New(correctionMgr->reconstructorMgr(), seedMgr.ptr(), timeCom, myCpu);
  fshMgr->tailManagerIs(fstMgr.ptr());

  HalfSliceNetwork::Ptr network = new HalfSliceNetwork(
      vectorSize,
      mapping.ptr(),
      myCpu,
      seedMgr.ptr(),
      hsMgr.ptr(), fshMgr.ptr(), fstMgr.ptr(),
      readerMgr.ptr(), writerMgr.ptr(),
      surrogate.ptr(), // HACK
      schedule.ptr(),
      seedInitializer.ptr(),
      (correctionStrategy == HalfSliceCorrectionNetworkImpl::NON_HOMOGENEOUS) ? seedInitializer.ptr() : NULL);

  double tac = getTime();
  log() << "Total Init Time = " << (tac - tic) / 1000.0 << " s\n";
  
  activityMgr->targetPhaseIs(lastIteration, PhaseRank(0));

  double toc = getTime();

  log() << "Total Solve Time = " << (toc - tic) / 1000.0 << " s\n";
  log() << "End LinPita\n";
}

} // end namespace Pita

// Entrypoint
Pita::LinearDriver::Ptr linearPitaDriverNew(SingleDomainDynamic<double> * pbDesc) {
  return Pita::LinearDriverImpl::New(pbDesc);
}
