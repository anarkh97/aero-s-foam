#include "Fwk.h"
#include "Types.h"

#include "ReducedLinearDriverImpl.h"

#include "DynamState.h"
#include <Problems.d/DynamDescr.h>

#include "LinearGenAlphaIntegrator.h"
#include "HomogeneousGenAlphaIntegrator.h"
#include "IntegratorPropagator.h"
#include "IntegratorSeedInitializer.h"

#include "StaticSliceStrategy.h"
#include "SliceMapping.h"

#include "LinearHalfSliceSchedule.h"
#include "HalfSliceCorrectionNetworkImpl.h"

#include "Seed.h"

#include "AffinePostProcessor.h"
#include "LinearFineIntegratorManager.h"
#include "HalfSlicePropagatorManager.h"
#include "HalfTimeSliceImpl.h"
#include "ReducedFullTimeSliceImpl.h"

#include <Comm.d/Communicator.h>
extern Communicator * structCom;

#include <Driver.d/Domain.h>
extern Domain * domain;

namespace Pita {

using namespace Hts;

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

  /* Available cpus */ 
  timeCom_ = structCom;
  remoteCoarse_ = false; //domain->solInfo().remoteCoarse && (timeCom_->numCPUs() > 1);
  numCpus_ = CpuCount(timeCom_->numCPUs() - (remoteCoarse_ ? 1 : 0));
  myCpu_ = CpuRank(timeCom_->myID());
 
  /* Other paramters */ 
  maxActive_ = HalfSliceCount(domain->solInfo().numTSperCycleperCPU);
  lastIteration_ = IterationRank(domain->solInfo().kiter);
  projectorTolerance_ = domain->solInfo().pitaProjTol;
  noForce_ = domain->solInfo().NoForcePita;
  
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
  /* Summarize problem and parameters */ 
  log() << "vectorSize = " << vectorSize_ << "\n";
  log() << "Slices = " << numSlices_ << ", MaxActive = " << maxActive_ << ", Cpus = " << numCpus_ << "\n";
  log() << "Num iter = " << lastIteration_ << "\n"; 
  log() << "dt = " << fineTimeStep_ << ", J/2 = " << halfSliceRatio_ << ", Dt = J*dt = " << coarseTimeStep_ << ", tf = " << finalTime_ << "\n";

  SliceMapping::Ptr mapping = SliceMapping::New(fullTimeSlices_, numCpus_, maxActive_.value(), StaticSliceStrategy::New().ptr());
  
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
  std::vector<int> localFileId;
  for (SliceMapping::SliceIdIterator it = mapping->hostedSlice(myCpu_, HalfSliceRank(0), HalfSliceRank(0) + mapping->totalSlices()); it; ++it) {
    if (it->type() == Hs::FORWARD_HALF_SLICE || it->type() == Hs::BACKWARD_HALF_SLICE) {
      localFileId.push_back(it->rank().value());
    }
  }
  std::sort(localFileId.begin(), localFileId.end());
  localFileId.erase(std::unique(localFileId.begin(), localFileId.end()), localFileId.end());
  
  AffinePostProcessor::Ptr pitaPostProcessor = AffinePostProcessor::New(geoSource, localFileId.size(), &localFileId[0], probDesc_->getPostProcessor());
  typedef PostProcessing::IntegratorReactorImpl<AffinePostProcessor> LinearIntegratorReactor;
  PostProcessing::Manager::Ptr postProcessingMgr = PostProcessing::Manager::New(LinearIntegratorReactor::Builder::New(pitaPostProcessor.ptr()).ptr());

  /* Seeds */ 
  Seed::Manager::Ptr seedMgr = Seed::Manager::New();
  ReducedSeed::Manager::Ptr reducedSeedMgr = ReducedSeed::Manager::New();

  /* Correction */
  HalfSliceSchedule::Ptr schedule = LinearHalfSliceSchedule::New(numSlices_);
  HalfSliceCorrectionNetworkImpl::Ptr correctionMgr =
    HalfSliceCorrectionNetworkImpl::New(
      vectorSize_,
      timeCom_,
      myCpu_,
      schedule.ptr(),
      mapping.ptr(),
      dynamOps.ptr(),
      HalfSliceCorrectionNetworkImpl::HOMOGENEOUS,
      projectorTolerance_); 

  /* Fine integrators */
  FineIntegratorManager::Ptr fineIntegratorMgr = LinearFineIntegratorManager<AffineGenAlphaIntegrator>::New(dopsManager.ptr(), fineIntegrationParam);

  HalfSlicePropagatorManager::Ptr propagatorMgr =
    new HalfSlicePropagatorManager(
      correctionMgr->collector(),
      fineIntegratorMgr.ptr(), 
      postProcessingMgr.ptr(),
      halfSliceRatio_,
      initialTime_);

  /* Slices */
  HalfTimeSliceImpl::Manager::Ptr hsMgr = HalfTimeSliceImpl::Manager::New(propagatorMgr.ptr());
  
  ReducedFullTimeSliceHeadImpl::Manager::Ptr fshMgr = ReducedFullTimeSliceHeadImpl::Manager::New(
      correctionMgr->updateProjectorMgr(),
      coarsePropagator.ptr(),
      seedMgr.ptr(),
      reducedSeedMgr.ptr(),
      correctionMgr->updatedSeedAssemblerMgr(),
      timeCom_,
      myCpu_);

  ReducedFullTimeSliceTailImpl::Manager::Ptr fstMgr = ReducedFullTimeSliceTailImpl::Manager::New(
      correctionMgr->updatedSeedAssemblerMgr(),
      seedMgr.ptr(),
      timeCom_,
      myCpu_);


  /* Instantiate the stuff */

  for (HalfSliceRank rank = HalfSliceRank(0); rank <= HalfSliceRank(0) + numSlices_; rank = rank + HalfSliceCount(1)) {
    Seed::Ptr mainSeed = seedMgr->instanceNew(String("MS_") + toString(rank)); 
    Seed::Ptr leftSeed = seedMgr->instanceNew(String("LPS_") + toString(rank)); 
    Seed::Ptr rightSeed = seedMgr->instanceNew(String("RPS_") + toString(rank));

    ReducedSeed::Ptr correctionComp = reducedSeedMgr->instanceNew(String("CS_") + toString(rank)); 
  }

  for (HalfSliceRank rank = HalfSliceRank(0); rank < HalfSliceRank(0) + numSlices_; rank = rank + HalfSliceCount(1)) {
    HalfTimeSliceImpl::Ptr forwardSlice = hsMgr->instanceNew(HalfSliceId(rank, HalfTimeSlice::FORWARD));
    forwardSlice->seedIs(seedMgr->instance(String("MS_") + toString(rank)));
    forwardSlice->propagatedSeedIs(seedMgr->instance(String("LPS_") + toString(rank + HalfSliceCount(1))));

    HalfTimeSliceImpl::Ptr backwardSlice = hsMgr->instanceNew(HalfSliceId(rank, HalfTimeSlice::BACKWARD));
    forwardSlice->seedIs(seedMgr->instance(String("MS_") + toString(rank + HalfSliceCount(1))));
    forwardSlice->propagatedSeedIs(seedMgr->instance(String("RPS_") + toString(rank)));
  }

  for (HalfSliceRank rank = HalfSliceRank(0); rank < HalfSliceRank(0) + numSlices_; rank = rank + HalfSliceCount(1)) {
    ReducedFullTimeSliceHeadImpl::Ptr headSlice = fshMgr->instanceNew(rank);
    headSlice->leftPropagatedSeedIs(seedMgr->instance(String("LPS_") + toString(rank)));
    headSlice->rightPropagatedSeedIs(seedMgr->instance(String("RPS_") + toString(rank)));
    headSlice->updatedSeedIs(seedMgr->instance(String("MS_") + toString(rank)));
    headSlice->correctionIs(reducedSeedMgr->instance(String("CS_") + toString(rank)));

    ReducedFullTimeSliceTailImpl::Ptr tailSlice = fstMgr->instanceNew(rank);
    tailSlice->nextLeftPropagatedSeedIs(seedMgr->instance(String("LPS_") + toString(rank + HalfSliceCount(1))));
    tailSlice->nextUpdatedSeedIs(seedMgr->instance(String("MS_") + toString(rank + HalfSliceCount(1))));
    headSlice->correctionIs(reducedSeedMgr->instance(String("CS_") + toString(rank + HalfSliceCount(1))));
  }

  /* End */ 
  log() << "End Reduced HalfSlice Linear Pita\n";
}

} /* end namespace Pita */

/* Entrypoint */
/*Pita::LinearDriver::Ptr
linearPitaDriverNew(SingleDomainDynamic<double> * pbDesc) {
  return Pita::ReducedLinearDriverImpl::New(pbDesc);
}*/
