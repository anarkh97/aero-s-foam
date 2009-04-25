#include <OOPita.d/NlDynamPostProcessor.h>
#include <OOPita.d/PitaNonLinDynam.h>

namespace Pita {

// NlDynamPostProcessor implementation
  
NlDynamPostProcessor::NlDynamPostProcessor(PitaNonLinDynamic * pbDesc) :
  probDesc_(pbDesc),
  dummy_(pbDesc ? pbDesc->solVecInfo() : 0)
{}

NlDynamPostProcessor::~NlDynamPostProcessor() {
  probDesc_->closeOutputFiles();
}

void
NlDynamPostProcessor::statusIs(DynamPostProcessor::Status s) {
  if (status() == DynamPostProcessor::OPEN)
    probDesc_->closeOutputFiles();
  if (s == DynamPostProcessor::OPEN) {
    probDesc_->openOutputFiles(sliceRank().value());
  }
  setStatus(s);
}

void
NlDynamPostProcessor::sliceRankIs(SliceRank r) {
  setSliceRank(r);
  statusIs(DynamPostProcessor::OPEN);
}

void
NlDynamPostProcessor::lastNlOutputIs(Seconds time, TimeStepCount step, const DynamState & state, const GeomState * geom, const GenVector<double> & force) {
  if (status() == DynamPostProcessor::CLOSED)
    statusIs(DynamPostProcessor::OPEN);
  probDesc_->pitaDynamOutput(sliceRank().value(), const_cast<GeomState*>(geom), const_cast<Vector&>(state.velocity()),
      dummy_, time.value(), step.value(), const_cast<Vector&>(force), dummy_);
}

// IntegratorReactor implementation

NlDynamPostProcessor::IntegratorReactor::IntegratorReactor(NlDynamTimeIntegrator * notifier, NlDynamPostProcessor * parent) :
  DynamPostProcessor::IntegratorReactor(notifier),
  nlNotifier_(notifier),
  parent_(parent)
{}

void
NlDynamPostProcessor::IntegratorReactor::onCurrentCondition() {
  parent()->lastNlOutputIs(nlNotifier()->currentTime(), nlNotifier()->timeStepCount(), nlNotifier()->currentState(), nlNotifier()->geomState(), nlNotifier()->externalForce());
}

void
NlDynamPostProcessor::IntegratorReactor::onInitialCondition() {
  this->onCurrentCondition();
}

NlDynamPostProcessor::IntegratorReactor *
NlDynamPostProcessor::getNewIntegratorReactor(const DynamTimeIntegrator * notifier) {
  return new IntegratorReactor(dynamic_cast<NlDynamTimeIntegrator *>(const_cast<DynamTimeIntegrator *>(notifier)), this);
}

NlDynamPostProcessor::IntegratorReactor::Ptr
NlDynamPostProcessor::integratorReactorNew(NlDynamTimeIntegrator * notifier) {
  return getNewIntegratorReactor(notifier);
}

  
} // end namespace Pita
