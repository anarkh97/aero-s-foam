#ifndef PITA_SEEDUPDATER_H
#define PITA_SEEDUPDATER_H

#include "Seed.h"
#include "NamedTask.h"

namespace Pita {

class SeedUpdater : public NamedTask {
public:
  EXPORT_PTRINTERFACE_TYPES(SeedUpdater);
  //class Manager;
  typedef Fwk::GenManagerInterface<SeedUpdater *, String> Manager;

  virtual void iterationIs(IterationRank ir);

  virtual size_t vectorSize() const;

  /* Sources */
  const Seed * propagatedSeed() const { return propagatedSeed_.ptr(); }
  Seed * correction() const { return correction_.ptr(); } // Not guaranteed to be immutable
  
  virtual void propagatedSeedIs(const Seed * ps) { setPropagatedSeed(ps); }
  virtual void correctionIs(Seed * c) { setCorrection(c); } 

  /* Targets */
  Seed * updatedSeed() const { return updatedSeed_.ptr(); }
 
  virtual void updatedSeedIs(Seed * us) { setUpdatedSeed(us); }

protected:
  explicit SeedUpdater(const String & name) :
    NamedTask(name)
  {}

  void updateSeed();

  void setPropagatedSeed(const Seed * ps) { propagatedSeed_ = ps; }
  void setCorrection(Seed * c) { correction_ = c; }
  void setUpdatedSeed(Seed * us) { updatedSeed_ = us; }

private:
  Seed::PtrConst propagatedSeed_;
  Seed::Ptr correction_;
  Seed::Ptr updatedSeed_;

  DISALLOW_COPY_AND_ASSIGN(SeedUpdater);
};

} /* end namespace Pita */

#endif /* PITA_SEEDUPDATER_H */
