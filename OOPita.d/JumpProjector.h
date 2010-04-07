#ifndef PITA_JUMPPROJECTOR_H
#define PITA_JUMPPROJECTOR_H

#include "Fwk.h"
#include "Seed.h"
#include "NamedTask.h"

namespace Pita {

class JumpProjector : public NamedTask {
public:
  EXPORT_PTRINTERFACE_TYPES(JumpProjector);
  typedef Fwk::GenManagerInterface<JumpProjector *, String> Manager;

  virtual size_t reducedBasisSize() const = 0;

  /* Sources */
  const Seed * predictedSeed() const { return predictedSeed_.ptr(); }
  const Seed * actualSeed() const { return actualSeed_.ptr(); }

  virtual void predictedSeedIs(const Seed * ps) { setPredictedSeed(ps); }
  virtual void actualSeedIs(const Seed * as) { setActualSeed(as); }

  /* Results */
  const Seed * seedJump() const { return seedJump_.ptr(); }
  const ReducedSeed * reducedSeedJump() const { return reducedSeedJump_.ptr(); }

  Seed * seedJump() { return seedJump_.ptr(); }
  ReducedSeed * reducedSeedJump() { return reducedSeedJump_.ptr(); }

  virtual void reducedSeedJumpIs(ReducedSeed * rsj) { setReducedSeedJump(rsj); }
  virtual void seedJumpIs(Seed * sj) { setSeedJump(sj); }

protected:
  explicit JumpProjector(const String & name) :
    NamedTask(name)
  {}

  void setPredictedSeed(const Seed * ps) { predictedSeed_ = ps; }
  void setActualSeed(const Seed * as) { actualSeed_ = as; }
  void setSeedJump(Seed * sj) { seedJump_ = sj; }
  void setReducedSeedJump(ReducedSeed * rsj) { reducedSeedJump_ = rsj; }
    
private:
  Seed::PtrConst predictedSeed_;
  Seed::PtrConst actualSeed_;
  Seed::Ptr seedJump_;
  ReducedSeed::Ptr reducedSeedJump_;

  DISALLOW_COPY_AND_ASSIGN(JumpProjector);
};

} /* end namespace Pita */

#endif /* PITA_JUMPPROJECTOR_H */
