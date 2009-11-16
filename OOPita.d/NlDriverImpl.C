#include "NlDriverImpl.h"
#include "NlTimeSlice.h"
#include "RemoteTimeSlice.h"
#include "CommManagerImpl.h"
#include "NlSeedInitializerImpl.h"
#include "NlDynamTimeIntegrator.h"
#include "NlDynamPostProcessor.h"
#include "LinearizedPropagator.h"
#include "NlDynamOps.h"
#include "ProjectorPropagator.h"
#include "IntegratorPropagator.h"
#include "Activity.h"
#include "Log.h"

extern Communicator * structCom;

/*Pita::NlDriver::Ptr nlPitaDriverNew(Pita::PitaNonLinDynamic * problemDescriptor) {
  return new Pita::NlDriverImpl(problemDescriptor) ;
}*/

namespace Pita {

NlDriverImpl::NlDriverImpl(PitaNonLinDynamic * pbDesc) :
  NlDriver(pbDesc),
  timeCommunicator_(structCom)
{}

void
NlDriverImpl::solve() {
  probDesc()->preProcess();

  Activity::Manager::Ptr activityMgr = activityManagerInstance(); 
  
  TimeSliceMapping::Ptr sliceMapping_(TimeSliceMapping::New(SliceCount(probDesc()->getNumTS()), SliceCount(probDesc()->getNumTSonCPU()), CpuCount(timeCommunicator_->numCPUs())));
  
  duplicateFiles(sliceMapping_.ptr());
 
  CommManager::Ptr timeCommManager_(CommManagerImpl::New(sliceMapping_.ptr(), structCom));
  
  NlDynamTimeIntegrator::Ptr integrator_ = NlDynamTimeIntegrator::New(probDesc());
  NlDynamPostProcessor::Ptr postProcessor_ = NlDynamPostProcessor::New(probDesc());
  DynamPostProcessor::IntegratorReactor::Ptr ppReactor_ = postProcessor_->integratorReactorNew(integrator_.ptr());
  
  NlTimeSlice::Manager::Ptr localSliceMgr_ = NlTimeSlice::Manager::New(integrator_.ptr(), timeCommManager_.ptr());
  RemoteTimeSlice::Manager::Ptr remoteSliceMgr_ = RemoteTimeSlice::Manager::New(timeCommManager_.ptr(), vectorSize());
  
  NlDynamTimeIntegrator::Ptr coarseIntegrator_ = NlDynamTimeIntegrator::New(probDesc());
  coarseIntegrator_->timeStepSizeIs(Seconds(probDesc()->getCoarseDt()));
  
  SeedInitializer::Ptr seedInitializer_ = NlSeedInitializerImpl::New(probDesc(), coarseIntegrator_.ptr(), TimeStepCount(1));

  TimeSliceNetwork::Ptr sliceNetwork_ = TimeSliceNetwork::New( myCpuRank(), sliceMapping_.ptr(), localSliceMgr_.ptr(), remoteSliceMgr_.ptr(), seedInitializer_.ptr() );

  if (probDesc()->getKiter()) {
    IterationRank lastIteration(probDesc()->getKiter() - 1);
    sliceNetwork_->lastIterationIs(lastIteration);
    activityMgr->targetPhaseIs(lastIteration, PhaseRank::fineGrid()); 
  }

  log() << "End NlPita\n";
}

void
NlDriverImpl::duplicateFiles(const TimeSliceMapping * sliceMapping) {
  int localSlices = sliceMapping->slicesOnCpu(myCpuRank()).value();
  std::vector<int> ts;
  ts.reserve(localSlices);
  for (TimeSliceMapping::SliceIteratorConst s = sliceMapping->slices(myCpuRank()); s; ++s) {
    ts.push_back((*s).value()); 
  }
  geoSource->duplicateFilesForPita(localSlices, &ts[0]);
}

} // end namespace Pita
