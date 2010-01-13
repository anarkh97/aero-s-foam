#ifndef PITA_UPDATEDSEEDASSEMBLER_H
#define PITA_UPDATEDSEEDASSEMBLER_H

#include "Fwk.h"
#include "Types.h"
#include "Seed.h"
#include "NamedTask.h"

namespace Pita {

class UpdatedSeedAssembler : public NamedTask {
public:
  EXPORT_PTRINTERFACE_TYPES(UpdatedSeedAssembler);
  typedef Fwk::GenManagerInterface<UpdatedSeedAssembler *, String> Manager;

  virtual size_t reducedBasisSize() const = 0;
  virtual size_t vectorSize() const = 0;

  /* Targets */
  Seed * correction() const { return correction_.ptr(); }
  Seed * updatedSeed() const { return updatedSeed_.ptr(); }
 
  virtual void correctionIs(Seed * c) = 0; 
  virtual void updatedSeedIs(Seed * us) = 0;

  /* Sources */
  const Seed * propagatedSeed() const { return propagatedSeed_.ptr(); }
  const ReducedSeed * correctionComponents() const { return correctionComponents_.ptr(); }
  
  virtual void propagatedSeedIs(const Seed * ps) = 0;
  virtual void correctionComponentsIs(const ReducedSeed * cc) = 0;

protected:
  explicit UpdatedSeedAssembler(const String & name) :
    NamedTask(name)
  {}

  void setCorrection(Seed * c) { correction_ = c; }
  void setUpdatedSeed(Seed * us) { updatedSeed_ = us; }
  void setPropagatedSeed(const Seed * ps) { propagatedSeed_ = ps; }
  void setCorrectionComponents(const ReducedSeed * cc) { correctionComponents_ = cc; }

private:
  Seed::Ptr correction_;
  Seed::Ptr updatedSeed_;
  Seed::PtrConst propagatedSeed_;
  ReducedSeed::PtrConst correctionComponents_; 

  DISALLOW_COPY_AND_ASSIGN(UpdatedSeedAssembler);
};

} /* end namespace Pita */

#endif /* PITA_UPDATEDSEEDASSEMBLER_H */
