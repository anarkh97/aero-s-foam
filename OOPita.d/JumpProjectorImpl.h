#ifndef PITA_JUMPPROJECTORIMPL_H
#define PITA_JUMPPROJECTORIMPL_H

#include "JumpProjector.h"

#include "JumpBuilder.h"
#include "DynamStateBasis.h"

namespace Pita {

class JumpProjectorImpl : public JumpProjector {
public:
  EXPORT_PTRINTERFACE_TYPES(JumpProjectorImpl);
  class Manager;

  virtual size_t reducedBasisSize() const;

  virtual void predictedSeedIs(const Seed * ps);
  virtual void actualSeedIs(const Seed * as);
  virtual void seedJumpIs(Seed * j);
  virtual void iterationIs(IterationRank i);

  /* Reduced basis and operators */
  const DynamStateBasis * reducedBasis() const;
  const JumpBuilder * jumpBuilder() const { return jumpBuilder_.ptr(); }

protected:
  JumpProjectorImpl(const String & name, const Manager * manager, JumpBuilder * jumpBuilder);

  friend class Manager;

private:
  const Manager * manager_;
  JumpBuilder::Ptr jumpBuilder_;
};

/* JumpProjectorImpl::Manager definition */

class JumpProjectorImpl::Manager : public JumpProjector::Manager, private Fwk::GenManagerImpl<JumpProjectorImpl, String> {
public:
  EXPORT_PTRINTERFACE_TYPES(Manager);

  // Overriden members
  virtual JumpProjectorImpl * instance(const String & key) const { return Impl::instance(key); }
  virtual size_t instanceCount() const { return Impl::instanceCount(); }
  virtual JumpProjectorImpl * instanceNew(const String & key) { return Impl::instanceNew(key); } 
  virtual void instanceDel(const String & key) { Impl::instanceDel(key); }

  // Added members
  const DynamStateBasis * defaultReducedBasis() const { return defaultReducedBasis_.ptr(); }

  void defaultReducedBasisIs(const DynamStateBasis * drb) { defaultReducedBasis_ = drb; }

  static Ptr New(const DynamStateBasis * defaultReducedBasis) {
    return new Manager(defaultReducedBasis);
  }

protected:
  explicit Manager(const DynamStateBasis * drb);

  virtual JumpProjectorImpl * createNewInstance(const String & key);

private:
  typedef Fwk::GenManagerImpl<JumpProjectorImpl, String> Impl;
  
  DynamStateBasis::PtrConst defaultReducedBasis_;
};

inline
const DynamStateBasis *
JumpProjectorImpl::reducedBasis() const {
  return manager_->defaultReducedBasis();
}

} /* end namespace Pita */

#endif /* PITA_JUMPPROJECTORIMPL_H */
