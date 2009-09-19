#include "UpdatedSeedAssemblerImpl.h"

#include <memory>

namespace Pita {

/* UpdatedSeedAssemblerImpl */

UpdatedSeedAssemblerImpl::UpdatedSeedAssemblerImpl(const String & name, const Manager * manager) :
  UpdatedSeedAssembler(name),
  manager_(manager)
  //schedulingReactor_(NULL),
  //updateReactor_(NULL)
{}

/*void
UpdatedSeedAssemblerImpl::assemblyPhaseIs(PhaseRank ap) {
  //updateReactor_->notifier()->phaseIs(ap);
  setAssemblyPhase(ap);
}*/

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
  //schedulingReactor_->notifierIs(cc);
  setCorrectionComponents(cc);
}

void
UpdatedSeedAssemblerImpl::iterationIs(IterationRank ir) {
  if (!propagatedSeed() || !updatedSeed()) {
    return;
  }

  //log() << "Reading state from " << correctionComponents()->name() << "\n";

  const Vector & components = correctionComponents()->state();

  DynamState result(propagatedSeed()->state());
  int rbs = static_cast<int>(correctionBasis()->stateCount());
  //log() << "reducedBasisSize/correctionComponentSize == " << rbs << "/" << components.size() << "\n";
  for (int i = 0; i < rbs; ++i) {
    if (components[i] != 0.0) {
      result.displacement().linAdd(components[i], correctionBasis()->state(i).displacement());
      result.velocity().linAdd(components[i], correctionBasis()->state(i).velocity());
    }
  }

  updatedSeed()->stateIs(result);
  updatedSeed()->iterationIs(ir);
  updatedSeed()->statusIs(propagatedSeed()->status());
}

/* SchedulingReactor */

/*void
UpdatedSeedAssemblerImpl::SchedulingReactor::onState() {
  if (activity()->status() != Activity::scheduled) {
    activity()->iterationIs(activity()->currentIteration());
    activity()->statusIs(Activity::scheduled);
  }
}*/

/* UpdateReactor */

/*void
UpdatedSeedAssemblerImpl::UpdateReactor::onStatus() {
  if (!parent_->propagatedSeed() || !parent()->updatedSeed()) {
    return;
  }

  const DynamStateBasis * basis = parent_->correctionBasis();
  const Vector & components = parent_->correctionComponents()->state();

  DynamState result(parent_->propagatedSeed()->state());
  int rbs = static_cast<int>(basis->stateCount());
  for (int i = 0; i < rbs; ++i) {
    if (components[i] != 0.0) {
      result.displacement().linAdd(components[i], basis->state(i).displacement());
      result.velocity().linAdd(components[i], basis->state(i).velocity());
    }
  }

  parent()->updatedSeed()->stateIs(result);
}*/

/* Manager */

UpdatedSeedAssemblerImpl::Manager::Manager(const DynamStateBasis * dcb) :
  defaultCorrectionBasis_(dcb)
{}

UpdatedSeedAssemblerImpl *
UpdatedSeedAssemblerImpl::Manager::createNewInstance(const String & key) {
  //String activityName = String("USA_") + key;
  //Activity::Ptr activity = activityManagerInstance()->activityNew(activityName);

  String taskName = String("UpdatedSeedAssembler ") + key;
  std::auto_ptr<UpdatedSeedAssemblerImpl> newInstance(new UpdatedSeedAssemblerImpl(taskName, this));
  
  //newInstance->schedulingReactor_ = new SchedulingReactor(NULL, activity.ptr());
  //newInstance->updateReactor_ = new UpdateReactor(activity.ptr(), newInstance.get());

  return newInstance.release();
}

} /* end namespace Pita */
