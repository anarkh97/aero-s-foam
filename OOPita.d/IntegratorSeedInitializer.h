#ifndef PITA_INTEGRATORSEEDINITIALIZER_H
#define PITA_INTEGRATORSEEDINITIALIZER_H

#include "Fwk.h"
#include "Types.h"

#include "SeedInitializer.h"

#include "DynamTimeIntegrator.h"
#include "DynamState.h"
#include "DynamStatePlainBasis.h"

namespace Pita {

class IntegratorSeedInitializer : public SeedInitializer {
public:
  EXPORT_PTRINTERFACE_TYPES(IntegratorSeedInitializer);

  virtual DynamState initialSeed(SliceRank rank) const; // Overriden
  
  DynamTimeIntegrator * integrator() const { return integrator_.ptr(); }
  TimeStepCount timeStepsBetweenSeeds() const { return timeStepsBetweenSeeds_; }

  static Ptr New(DynamTimeIntegrator * integrator,
                 TimeStepCount timeStepsBetweenSeeds) {
    return new IntegratorSeedInitializer(integrator, timeStepsBetweenSeeds);
  }

protected:
  IntegratorSeedInitializer(DynamTimeIntegrator * i, TimeStepCount tsbs);

private:
  DynamTimeIntegrator::Ptr integrator_;
  TimeStepCount timeStepsBetweenSeeds_;

  DynamStatePlainBasis::Ptr seed_;
};

} // end namespace Pita

#endif /* PITA_INTEGRATORSEEDINITIALIZER_H */
