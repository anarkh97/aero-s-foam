#ifndef PITA_NLSEEDINITIALIZERIMPL_H
#define PITA_NLSEEDINITIALIZERIMPL_H

#include "Fwk.h"
#include "Types.h"

#include "SeedInitializer.h"

#include "DynamTimeIntegrator.h"
#include "DynamState.h"
#include "DynamStatePlainBasis.h"

namespace Pita {

class PitaNonLinDynamic;
  
class NlSeedInitializerImpl : public SeedInitializer {
public:
  EXPORT_PTRINTERFACE_TYPES(NlSeedInitializerImpl);

  virtual DynamState initialSeed(SliceRank rank) const; // Overriden
  
  DynamTimeIntegrator * integrator() const { return integrator_.ptr(); }
  TimeStepCount timeStepsBetweenSeeds() const { return timeStepsBetweenSeeds_; }

  static NlSeedInitializerImpl::Ptr New(PitaNonLinDynamic * pbDesc,
                                        DynamTimeIntegrator * integrator,
                                        TimeStepCount timeStepsBetweenSeeds) {
    return new NlSeedInitializerImpl(pbDesc, integrator, timeStepsBetweenSeeds);
  }

protected:
  NlSeedInitializerImpl(PitaNonLinDynamic * pd,
                        DynamTimeIntegrator * ti,
                        TimeStepCount tsbs);

private:
  PitaNonLinDynamic * probDesc_;

  SliceRank lastFileSlice_;

  TimeStepCount timeStepsBetweenSeeds_;

  mutable DynamTimeIntegrator::Ptr integrator_;
  mutable DynamStatePlainBasis::Ptr state_;
};

} // end namespace Pita

#endif /* PITA_SEEDINITIALIZERIMPL_H */
