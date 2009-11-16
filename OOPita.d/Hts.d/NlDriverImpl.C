#include "NlDriverImpl.h"

#include "../NlDynamTimeIntegrator.h"
#include "../IntegratorPropagator.h"

#include "../NlPostProcessor.h"
#include "../PostProcessingManager.h"

#include "NlTaskManager.h"

namespace Pita { namespace Hts {

NlDriverImpl::NlDriverImpl(PitaNonLinDynamic * probDesc,
                           GeoSource * geoSource,
                           Domain * domain,
                           SolverInfo * solverInfo,
                           Communicator * baseComm) :
  NlDriver(probDesc),
  geoSource_(geoSource),
  domain_(domain),
  solverInfo_(solverInfo),
  baseComm_(baseComm)
{}

void
NlDriverImpl::solve() {
  log() << "Begin Reversible NonLinear Pita\n";
  double tic = getTime(); // Total time

  preprocess();

  /* Summarize problem and parameters */
  log() << "\n"; 
  log() << "Slices = " << mapping_->totalSlices() << ", MaxActive = " << mapping_->maxWorkload() << ", Cpus = " << mapping_->availableCpus() << "\n";
  log() << "Iteration count = " << lastIteration_ << "\n"; 
  log() << "dt = " << fineTimeStep_ << ", J/2 = " << halfSliceRatio_ << ", Dt = J*dt = " << coarseTimeStep_ << ", Tf = Slices*(J/2)*dt = " << finalTime_ << "\n";
  log() << "VectorSize = " << vectorSize_ << " dofs\n";
  log() << "Projector tol = " << projectorTolerance_ << "\n";

  // Problem parameters for the sequential-in-time
  DynamState initialState(vectorSize_);
  probDesc()->getInitState(initialState);
  TimeStepCount timeStepCount = halfSliceRatio_ * mapping_->totalSlices().value();
  NlDynamTimeIntegrator::Ptr integrator = NlDynamTimeIntegrator::New(probDesc());
  
  // Propagators
  IntegratorPropagator::Ptr forwardPropagator = IntegratorPropagator::New(integrator.ptr());
  forwardPropagator->initialTimeIs(Seconds(0.0));
  forwardPropagator->timeStepSizeIs(fineTimeStep_);
  forwardPropagator->timeStepCountIs(timeStepCount);
  
  IntegratorPropagator::Ptr backwardPropagator = IntegratorPropagator::New(integrator.ptr());
  backwardPropagator->initialTimeIs(finalTime_);
  backwardPropagator->timeStepSizeIs(-fineTimeStep_);
  backwardPropagator->timeStepCountIs(timeStepCount);
 
  // Set-up post-processing 
  /*CpuRank myCpu(baseComm()->myID()); 
  std::vector<int> ts;
  for (SliceMapping::SliceIterator s = mapping_->hostedSlice(myCpu); s; ++s) {
    ts.push_back((*s).value()); 
  }*/
  const int outputFileCount = 2;
  int ts[outputFileCount] = { 0, 1 };
  NlPostProcessor::Ptr pitaPostProcessor = NlPostProcessor::New(geoSource(), outputFileCount, ts, probDesc());
  typedef PostProcessing::IntegratorReactorImpl<NlPostProcessor> NlIntegratorReactor;
  NlIntegratorReactor::Builder::Ptr ppBuilder = NlIntegratorReactor::Builder::New(pitaPostProcessor.ptr());
  PostProcessing::Manager::Ptr ppMgr = PostProcessing::Manager::New(ppBuilder.ptr());
  
  ppMgr->outputFileSetIs(forwardPropagator.ptr(), PostProcessor::FileSetId(0));
  ppMgr->outputFileSetIs(backwardPropagator.ptr(), PostProcessor::FileSetId(1));
 
  // Perform forward-in-time then backward-in-time
  forwardPropagator->initialStateIs(initialState);
  backwardPropagator->initialStateIs(forwardPropagator->finalState());
  
  double toc = getTime();
  log() << "\n" << "Total time = " << (toc - tic) / 1000.0 << " s\n";
  log() << "\n" << "End Reversible NonLinear Pita\n";
}

void
NlDriverImpl::preprocess() {
  double tic = getTime();
  
  probDesc()->preProcess();
  
  /* Space-domain */
  vectorSize_ = probDesc()->solVecInfo();

  /* Time-domain */
  fineTimeStep_ = Seconds(solverInfo()->dt);
  halfSliceRatio_ = TimeStepCount(solverInfo()->Jratio / 2);
  sliceRatio_ = TimeStepCount(halfSliceRatio_.value() * 2);
  coarseTimeStep_ = fineTimeStep_ * sliceRatio_.value(); 
  initialTime_ = Seconds(solverInfo()->initialTime);
  finalTime_ = Seconds(solverInfo()->tmax);
 
  HalfSliceCount numSlices(static_cast<int>(ceil((finalTime_.value() - initialTime_.value()) / (halfSliceRatio_.value() * fineTimeStep_.value()))));
  FullSliceCount fullTimeSlices((numSlices.value() / 2) + (numSlices.value() % 2));
  numSlices = HalfSliceCount(fullTimeSlices.value() * 2);
  finalTime_ = fineTimeStep_ * Seconds(numSlices.value() * halfSliceRatio_.value()); // To have a whole number of full time-slices 

  /* Load balancing */ 
  CpuCount numCpus(baseComm()->numCPUs());
  HalfSliceCount maxActive(solverInfo()->numTSperCycleperCPU);
  mapping_ = SliceMapping::New(fullTimeSlices, numCpus, maxActive);

  /* Other parameters */ 
  lastIteration_ = IterationRank(solverInfo()->kiter);
  projectorTolerance_ = solverInfo()->pitaProjTol;
 
  double toc = getTime();
  log() << "\n";
  log() << "Total preprocessing time = " << (toc - tic) / 1000.0 << " s\n";
}

} /* end namespace Hts */ } /* end namespace Pita */

extern GeoSource * geoSource;
extern Domain * domain;
extern Communicator * structCom;

Pita::NlDriver::Ptr
nlPitaDriverNew(Pita::PitaNonLinDynamic * problemDescriptor) {
  return Pita::Hts::NlDriverImpl::New(problemDescriptor, geoSource, domain, &domain->solInfo(), structCom);
}

