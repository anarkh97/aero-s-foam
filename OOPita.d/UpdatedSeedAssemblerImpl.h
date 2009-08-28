#ifndef PITA_UPDATEDSEEDASSEMBLERIMPL_H
#define PITA_UPDATEDSEEDASSEMBLERIMPL_H

#include "UpdatedSeedAssembler.h"

#include "DynamStateBasis.h"
#include "Activity.h"

namespace Pita {

class UpdatedSeedAssemblerImpl : public UpdatedSeedAssembler {
public:
  EXPORT_PTRINTERFACE_TYPES(UpdatedSeedAssemblerImpl);
  class Manager;

  virtual size_t reducedBasisSize() const { return correctionBasis()->stateCount(); }

  using UpdatedSeedAssembler::updatedSeed;
  Seed * updatedSeed() { return const_cast<Seed *>(const_cast<const UpdatedSeedAssemblerImpl *>(this)->updatedSeed()); }

  /* Overriden members */
  virtual void assemblyPhaseIs(PhaseRank ap);
  virtual void updatedSeedIs(Seed * us);
  virtual void propagatedSeedIs(const Seed * ps);
  virtual void correctionComponentsIs(const ReducedSeed * cc);

  /* Added members */
  const Manager * manager() const { return manager_; }
  const DynamStateBasis * correctionBasis() const;

  virtual void doAssembly();

protected:
  friend class Manager;
  class SchedulingReactor;
  class UpdateReactor;

  explicit UpdatedSeedAssemblerImpl(const Manager * manager);

private:
  const Manager * manager_; // Back pointer
  Fwk::Ptr<SchedulingReactor> schedulingReactor_;
  Fwk::Ptr<UpdateReactor> updateReactor_;
};


class UpdatedSeedAssemblerImpl::SchedulingReactor : public ReducedSeed::NotifieeConst {
public:
  EXPORT_PTRINTERFACE_TYPES(SchedulingReactor);

  virtual void onState(); // overriden

  const Activity * activity() const { return activity_.ptr(); }
  Activity * activity() { return activity_.ptr(); }
  void activityIs(Activity * a) { activity_ = a; }

  SchedulingReactor(const ReducedSeed * n, Activity * a) :
    ReducedSeed::NotifieeConst(n),
    activity_(a)
  {}

private:
  Activity::Ptr activity_;
};

class UpdatedSeedAssemblerImpl::UpdateReactor : public Activity::Notifiee {
public:
  EXPORT_PTRINTERFACE_TYPES(UpdateReactor);

  // overriden
  virtual void onStatus();

  const UpdatedSeedAssemblerImpl * parent() const { return parent_; }
  UpdatedSeedAssemblerImpl * parent() { return parent_; }

  UpdateReactor(Activity * notifier, UpdatedSeedAssemblerImpl * parent) :
    Activity::Notifiee(notifier),
    parent_(parent)
  {}

private:
  UpdatedSeedAssemblerImpl * parent_; // Back pointer
};

class UpdatedSeedAssemblerImpl::Manager : public UpdatedSeedAssembler::Manager,
                                          private Fwk::GenManagerImpl<UpdatedSeedAssemblerImpl, String> {
public:
  EXPORT_PTRINTERFACE_TYPES(Manager);

  /* Overiden members */
  virtual UpdatedSeedAssemblerImpl * instance(const String & key) const { return Impl::instance(key); }
  virtual size_t instanceCount() const { return Impl::instanceCount(); }
  virtual UpdatedSeedAssemblerImpl * instanceNew(const String & key) { return Impl::instanceNew(key); } 
  virtual void instanceDel(const String & key) { Impl::instanceDel(key); }

  /* Added members */
  const DynamStateBasis * defaultCorrectionBasis() const { return defaultCorrectionBasis_.ptr(); }
  void defaultCorrectionBasisIs(const DynamStateBasis * dcb) { defaultCorrectionBasis_ = dcb; }

  static Ptr New(const DynamStateBasis * defaultCorrectionBasis) {
    return new Manager(defaultCorrectionBasis);
  }

protected:
  explicit Manager(const DynamStateBasis * defaultCorrectionBasis);

private:
  typedef Fwk::GenManagerImpl<UpdatedSeedAssemblerImpl, String> Impl;
  
  // Overriden
  virtual UpdatedSeedAssemblerImpl * createNewInstance(const String & key);

  DynamStateBasis::PtrConst defaultCorrectionBasis_;
};

inline
const DynamStateBasis *
UpdatedSeedAssemblerImpl::correctionBasis() const {
  return manager_->defaultCorrectionBasis();
}

} /* end namespace Pita */

#endif /* PITA_UPDATEDSEEDASSEMBLERIMPL_H */
