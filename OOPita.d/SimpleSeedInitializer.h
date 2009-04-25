#ifndef PITA_HTS_SIMPLESEEDINITIALIZER_H
#define PITA_HTS_SIMPLESEEDINITIALIZER_H

#include "SeedInitializer.h"

namespace Pita { namespace Hts {
class SimpleSeedInitializer : public SeedInitializer {
public:
  EXPORT_PTRINTERFACE_TYPES(SimpleSeedInitializer);

  virtual DynamState initialSeed(SliceRank rank) const; // Overriden

  static Ptr New(const DynamState & firstInitialSeed) {
    return new SimpleSeedInitializer(firstInitialSeed);
  }

protected:
  explicit SimpleSeedInitializer(const DynamState & firstInitialSeed);

private:
  DynamState firstInitialSeed_;
};

} /* end namespace Hts */ } /* end namespace Pita */

#endif /* PITA_HTS_SIMPLESEEDINITIALIZER_H */
