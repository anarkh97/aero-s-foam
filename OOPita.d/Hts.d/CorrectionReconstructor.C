#include "CorrectionReconstructor.h"

namespace Pita { namespace Hts {

void
CorrectionReconstructor::iterationIs(IterationRank ir) {
  if (correction()->iteration() < correctionComponents()->iteration()) {
    if (correctionComponents()->status() != Seed::INACTIVE) {
      log() << "*** Assembling from projection on space of size " << reconstructor_->reducedBasisSize() << "\n";
      reconstructor_->reducedBasisComponentsIs(correctionComponents()->state());
      correction()->stateIs(reconstructor_->finalState());
    }
    correction()->statusIs(correctionComponents()->status());
    correction()->iterationIs(correctionComponents()->iteration());
  }
}

CorrectionReconstructor::Manager::Manager(OperatorManager * reconstructorManager) :
  impl_(Factory(reconstructorManager))
{}

CorrectionReconstructor *
CorrectionReconstructor::Manager::Factory::operator()(HalfSliceRank key) {
  DynamStateReconstructor * reconstructor = manager_->instanceNew(key);
  String taskName = String("CorrectionReconstructor ") + toString(key);
  return new CorrectionReconstructor(taskName, reconstructor);
}

} /* end namespace Hts */ } /* end namespace Pita */
