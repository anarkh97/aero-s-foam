#ifndef PITA_JUMPPROJECTOR_H
#define PITA_JUMPPROJECTOR_H

#include "Fwk.h"
#include "Seed.h"
#include "JumpBuilder.h"

namespace Pita {

class JumpProjector : public JumpBuilder {
public:
  EXPORT_PTRINTERFACE_TYPES(JumpProjector);
  typedef Fwk::GenManagerSubInterface<JumpProjector *, String, Fwk::GenManagerInterface<JumpBuilder *, String> > Manager;

  virtual size_t reducedBasisSize() const = 0;

  /* Results */
  const ReducedSeed * reducedSeedJump() const { return reducedSeedJump_.ptr(); }
  ReducedSeed * reducedSeedJump() { return reducedSeedJump_.ptr(); }

  virtual void reducedSeedJumpIs(ReducedSeed * rsj) { setReducedSeedJump(rsj); }

protected:
  explicit JumpProjector(const String & name) :
    JumpBuilder(name)
  {}

  void setReducedSeedJump(ReducedSeed * rsj) { reducedSeedJump_ = rsj; }
    
private:
  ReducedSeed::Ptr reducedSeedJump_;

  DISALLOW_COPY_AND_ASSIGN(JumpProjector);
};

} /* end namespace Pita */

#endif /* PITA_JUMPPROJECTOR_H */
