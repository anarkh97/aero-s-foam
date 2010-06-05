#ifndef PITA_SEEDERROREVALUATOR_H
#define PITA_SEEDERROREVALUATOR_H

#include "Fwk.h"
#include "Seed.h"
#include "DynamOps.h"

namespace Pita {

class SeedErrorEvaluator : public SharedState<DynamState>::NotifieeConst {
public:
  EXPORT_PTRINTERFACE_TYPES(SeedErrorEvaluator);

  class Manager;

  const DynamOps * metric() const { return metric_.ptr(); }
  const Seed * seedError() const { return notifier(); } 
  const Seed * referenceSeed() const { return referenceSeed_.ptr(); }

  void referenceSeedIs(const Seed * rs) { referenceSeed_ = rs; }

  // overriden
  virtual void onIteration();

protected:
  SeedErrorEvaluator(const DynamOps * metric, const Seed * seedError, const Seed * referenceSeed);

  friend class Manager;

private:
  DynamOps::PtrConst metric_;
  Seed::PtrConst referenceSeed_;
};

class SeedErrorEvaluator::Manager : public Fwk::GenManager<SeedErrorEvaluator, const Seed *> {
public:
  EXPORT_PTRINTERFACE_TYPES(Manager);

  const DynamOps * defaultMetric() const { return defaultMetric_.ptr(); }

  static Ptr New(const DynamOps * defaultMetric) {
    return new Manager(defaultMetric);
  }

protected:
  explicit Manager(const DynamOps * defaultMetric);

  virtual SeedErrorEvaluator * createNewInstance(const KeyType & key); // overriden

private:
  DynamOps::PtrConst defaultMetric_;
};

} /* end namespace Pita */

#endif /* PITA_SEEDERROREVALUATOR_H */
