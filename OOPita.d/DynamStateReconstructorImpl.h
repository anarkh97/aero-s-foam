#ifndef PITA_DYNAMSTATERECONSTRUCTORIMPL_H
#define PITA_DYNAMSTATERECONSTRUCTORIMPL_H

#include "DynamStateReconstructor.h"

#include "DynamStateBasis.h"

namespace Pita {

class DynamStateReconstructorImpl : public DynamStateReconstructor {
public:
  EXPORT_PTRINTERFACE_TYPES(DynamStateReconstructorImpl);

  // Overriden members
  virtual void reducedBasisComponentsIs(const Vector & c);

  const DynamStateBasis * reconstructionBasis() const { return reconstructionBasis_.ptr(); } 

  class Manager;

protected:
  explicit DynamStateReconstructorImpl(const DynamStateBasis * rb);

  void setReconstructionBasis(const DynamStateBasis * rb) { reconstructionBasis_ = rb; }

  void reset();

private:
  DynamStateBasis::PtrConst reconstructionBasis_;
};

class DynamStateReconstructorImpl::Manager : public DynamStateReconstructor::Manager,
                                             private Fwk::GenManagerImpl<DynamStateReconstructorImpl, Fwk::String> {
public:
  EXPORT_PTRINTERFACE_TYPES(Manager);

  // Overriden members
  virtual DynamStateReconstructorImpl * instance(const String & key) const { return Impl::instance(key); }
  virtual size_t instanceCount() const { return Impl::instanceCount(); }
  virtual DynamStateReconstructorImpl * instanceNew(const String & key) { return Impl::instanceNew(key); }
  virtual void instanceDel(const String & key) { Impl::instanceDel(key); }

  // Added members
  const DynamStateBasis * defaultReconstructionBasis() const { return defaultReconstructionBasis_.ptr(); }
  void defaultReconstructionBasisIs(const DynamStateBasis * drb);

  static Ptr New(const DynamStateBasis * drb) {
    return new Manager(drb);
  }

protected:
  explicit Manager(const DynamStateBasis * drb);

  virtual DynamStateReconstructorImpl * createNewInstance(const String & key);

  void setDefaultReconstructionBasis(const DynamStateBasis * drb) { defaultReconstructionBasis_ = drb; }

  void resetInstances();  

private:
  typedef Fwk::GenManagerImpl<DynamStateReconstructorImpl, Fwk::String> Impl;

  DynamStateBasis::PtrConst defaultReconstructionBasis_;
};

}

#endif /* PITA_DYNAMSTATERECONSTRUCTORIMPL_H */
