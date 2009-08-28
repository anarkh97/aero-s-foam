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
  
  const Seed * seed(SeedId id) const;
  void seedIs(SeedId id, const Seed * s);

  SeedStatus seedStatus(SeedId id) const;
  void seedStatusIs(SeedId id, SeedStatus status);

  const Seed * jumpSeed() const { return jumpSeed_.ptr(); }
  Seed * jumpSeed() { return jumpSeed_.ptr(); }
  void jumpSeedIs(Seed * j);

  static Ptr New() {
    return new JumpBuilder();
  }

protected:
  class SeedReactor;

  JumpBuilder();

private:
  SeedStatus seedStatus_[2];
  Fwk::Ptr<SeedReactor> seedReactor_[2];
  Seed::Ptr jumpSeed_;
};

} /* end namespace Pita */

#endif /* PITA_JUMPBUILDER_H */
