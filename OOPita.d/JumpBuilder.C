#include "JumpBuilder.h"

namespace Pita {

/* JumpBuilder::SeedReactor definition */

class JumpBuilder::SeedReactor : public Seed::NotifieeConst {
public:
  EXPORT_PTRINTERFACE_TYPES(SeedReactor);

  virtual void onState(); // overriden

  SeedReactor(const Seed * notifier, JumpBuilder * parent, SeedId id);

private:
  JumpBuilder * parent_;
  SeedId id_;
};

/* JumpBuilder::SeedReactor implementation */

JumpBuilder::SeedReactor::SeedReactor(
    const Seed * notifier,
    JumpBuilder * parent,
    SeedId id) :
  Seed::NotifieeConst(notifier),
  parent_(parent),
  id_(id)
{}

void
JumpBuilder::SeedReactor::onState() {
  parent_->seedStatusIs(id_, UPDATED);
}

/* JumpBuilder implementation */

JumpBuilder::JumpBuilder() :
  seedStatus_(),
  seedReactor_(),
  jump_()
{
  seedStatus_[LEFT] = seedStatus_[RIGHT] = MISSING;
  seedReactor_[LEFT] = new SeedReactor(NULL, this, LEFT);
  seedReactor_[RIGHT] = new SeedReactor(NULL, this, RIGHT);
}

const Seed *
JumpBuilder::propagatedSeed(JumpBuilder::SeedId id) const {
  return seedReactor_[id]->notifier();
}

void
JumpBuilder::propagatedSeedIs(JumpBuilder::SeedId id, const Seed * s) {
  if (propagatedSeed(id) != s) {
    seedReactor_[id]->notifierIs(s);
    SeedStatus newStatus = seedReactor_[id]->notifier() ? CURRENT : MISSING;
    seedStatusIs(id, newStatus);
  }
}

JumpBuilder::SeedStatus
JumpBuilder::seedStatus(JumpBuilder::SeedId id) const {
  return seedStatus_[id];
}

void
JumpBuilder::seedStatusIs(JumpBuilder::SeedId id, JumpBuilder::SeedStatus status) {
  if (seedStatus(id) != status) {
    seedStatus_[id] = status;
    if (seedStatus_[RIGHT] == UPDATED && seedStatus_[LEFT] == UPDATED) {
      jump_ = propagatedSeed(LEFT)->state() - propagatedSeed(RIGHT)->state();
      seedStatus_[RIGHT] = seedStatus_[LEFT] = CURRENT;
      if (notifiee_) {
        notifiee_->onJump();
      }
    }
  } 
}

} /* end namespace Pita */
