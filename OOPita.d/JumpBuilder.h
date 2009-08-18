#ifndef PITA_JUMPBUILDER_H
#define PITA_JUMPBUILDER_H

#include "Fwk.h"
#include "Seed.h"

namespace Pita {

class JumpBuilder : public Fwk::PtrInterface<JumpBuilder> {
public:
  EXPORT_PTRINTERFACE_TYPES(JumpBuilder);
  
  enum SeedId {
    LEFT = 0,
    RIGHT = 1
  };

  enum SeedStatus {
    MISSING = 0,
    CURRENT,
    UPDATED
  };
  
  const Seed * propagatedSeed(SeedId id) const;
  void propagatedSeedIs(SeedId id, const Seed * s);

  SeedStatus seedStatus(SeedId id) const;
  void seedStatusIs(SeedId id, SeedStatus status);

  const DynamState & jump() const { return jump_; }

  static Ptr New() {
    return new JumpBuilder();
  }

  class NotifieeConst : public Fwk::BaseNotifiee<const JumpBuilder, NotifieeConst> {
  public:
    EXPORT_PTRINTERFACE_TYPES(NotifieeConst);

    virtual void onJump() {}

  protected:
    explicit NotifieeConst(const JumpBuilder * notifier) :
      Fwk::BaseNotifiee<const JumpBuilder, NotifieeConst>(notifier)
    {}
  };

  NotifieeConst::Ptr lastNotifiee() const { return notifiee_; }
  void lastNotifieeIs(NotifieeConst * n) const { notifiee_ = n; }

protected:
  class SeedReactor;

  JumpBuilder();

private:
  SeedStatus seedStatus_[2];
  Fwk::Ptr<SeedReactor> seedReactor_[2];
  DynamState jump_;

  mutable NotifieeConst::Ptr notifiee_;
};

} /* end namespace Pita */

#endif /* PITA_JUMPBUILDER_H */
