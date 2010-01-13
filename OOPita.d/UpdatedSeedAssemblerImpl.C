#include "UpdatedSeedAssemblerImpl.h"

#include <memory>
#include <cassert>

namespace Pita {

/* UpdatedSeedAssemblerImpl */

UpdatedSeedAssemblerImpl::UpdatedSeedAssemblerImpl(const String & name, const Manager * manager) :
  UpdatedSeedAssembler(name),
  manager_(manager)
{}

void
UpdatedSeedAssemblerImpl::correctionIs(Seed * c) {
  setCorrection(c);
}

void
UpdatedSeedAssemblerImpl::updatedSeedIs(Seed * us) {
  setUpdatedSeed(us);
}

void
UpdatedSeedAssemblerImpl::propagatedSeedIs(const Seed * ps) {
  setPropagatedSeed(ps);
}

void
UpdatedSeedAssemblerImpl::correctionComponentsIs(const ReducedSeed * cc) {
  setCorrectionComponents(cc);
}

void
UpdatedSeedAssemblerImpl::iterationIs(IterationRank ir) {
  if (!propagatedSeed() || !updatedSeed()) {
    return;
  }

  if (correction()->iteration() < correctionComponents()->iteration()) {
  
    //log() << "Reading state from " << correctionComponents()->name() << "\n";
    const Vector & components = correctionComponents()->state();
    
    DynamState result(propagatedSeed()->state().vectorSize(), 0.0);
    int rbs = static_cast<int>(correctionBasis()->stateCount());
    //log() << "reducedBasisSize/correctionComponentSize == " << rbs << "/" << components.size() << "\n";
    for (int i = 0; i < rbs; ++i) {
      if (components[i] != 0.0) {
        result.linAdd(components[i], correctionBasis()->state(i));
      }
    }

    correction()->stateIs(result);
    correction()->statusIs(correctionComponents()->status());
    correction()->iterationIs(correctionComponents()->iteration());
  }

  assert(correction()->iteration() == ir);
  assert(propagatedSeed()->iteration() == ir);
  
  updatedSeed()->stateIs(correction()->state() + propagatedSeed()->state());
  updatedSeed()->statusIs(propagatedSeed()->status());
  updatedSeed()->iterationIs(propagatedSeed()->iteration());
}

/* Manager */

UpdatedSeedAssemblerImpl::Manager::Manager(const DynamStateBasis * dcb) :
  defaultCorrectionBasis_(dcb)
{}

UpdatedSeedAssemblerImpl *
UpdatedSeedAssemblerImpl::Manager::createNewInstance(const String & key) {
  String taskName = String("UpdatedSeedAssembler ") + key;
  std::auto_ptr<UpdatedSeedAssemblerImpl> newInstance(new UpdatedSeedAssemblerImpl(taskName, this));
  return newInstance.release();
}

} /* end namespace Pita */
