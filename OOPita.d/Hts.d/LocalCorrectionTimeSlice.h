#ifndef PITA_HTS_LOCALCORRECTIONTIMESLICE_H
#define PITA_HTS_LOCALCORRECTIONTIMESLICE_H

#include "CorrectionTimeSlice.h"

//#include "../JumpBuilder.h"
#include "../DynamPropagator.h"

namespace Pita { namespace Hts {

class LocalCorrectionTimeSlice : public CorrectionTimeSlice {
public:
  EXPORT_PTRINTERFACE_TYPES(LocalCorrectionTimeSlice);
  class Manager;

  virtual void predictedSeedIs(const Seed * ps); // overriden
  virtual void actualSeedIs(const Seed * as); // overriden
  virtual void jumpIs(Seed * j); // overriden

  virtual void iterationIs(IterationRank i); // overriden

protected:
  friend class Manager;
 
  LocalCorrectionTimeSlice(HalfSliceRank head, DynamPropagator * prop);
  
private:
  //JumpBuilder::Ptr jumpBuilder_;
  DynamPropagator::Ptr propagator_;  
};

class LocalCorrectionTimeSlice::Manager : public CorrectionTimeSlice::Manager, private Fwk::GenManagerImpl<LocalCorrectionTimeSlice, HalfSliceRank> {
public:
  EXPORT_PTRINTERFACE_TYPES(Manager);

  // Instance management
  LocalCorrectionTimeSlice * instance(const HalfSliceRank & key) const { return Impl::instance(key); }
  InstanceCount instanceCount() const { return Impl::instanceCount(); }
  LocalCorrectionTimeSlice * instanceNew(const HalfSliceRank & key) { return Impl::instanceNew(key); }
  void instanceDel(const HalfSliceRank & key) { Impl::instanceDel(key); }
  
  // Propagator is shared by all instances
  DynamPropagator * sharedPropagator() const { return sharedPropagator_.ptr(); }
  void sharedPropagatorIs(DynamPropagator * sp) { sharedPropagator_ = sp; }
  
  static Ptr New(DynamPropagator * sharedPropagator) {
    return new Manager(sharedPropagator);
  }
  
protected:
  explicit Manager(DynamPropagator * sp);
  
  virtual LocalCorrectionTimeSlice * createNewInstance(const HalfSliceRank & key); // overriden
  
private:
  typedef Fwk::GenManagerImpl<LocalCorrectionTimeSlice, HalfSliceRank> Impl;
  
  DynamPropagator::Ptr sharedPropagator_;

  DISALLOW_COPY_AND_ASSIGN(Manager);
};

} /* end namespace Hts */ } /* end namespace Pita */

#endif /* PITA_HTS_LOCALCORRECTIONTIMESLICE_H */
