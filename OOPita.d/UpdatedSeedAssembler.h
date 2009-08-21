#ifndef PITA_UPDATEDSEEDASSEMBLER_H
#define PITA_UPDATEDSEEDASSEMBLER_H

#include "Fwk.h"
#include "Types.h"
#include "Seed.h"

namespace Pita {

class UpdatedSeedAssembler : public Fwk::PtrInterface<UpdatedSeedAssembler> {
public:
  EXPORT_PTRINTERFACE_TYPES(UpdatedSeedAssembler);
  typedef Fwk::GenManagerInterface<UpdatedSeedAssembler *, String> Manager;

  /* Phase */
  PhaseRank assemblyPhase() const { return assemblyPhase_; }

  virtual void assemblyPhaseIs(PhaseRank ap) = 0;

  /* Client */
  const Seed * updatedSeed() const { return updatedSeed_.ptr(); }
  
  virtual void updatedSeedIs(Seed * us) = 0;

  /* Sources */
  const Seed * propagatedSeed() const { return propagatedSeed_.ptr(); }
  const ReducedSeed * correctionComponents() const { return correctionComponents_.ptr(); }
  
  virtual void propagatedSeedIs(const Seed * ps) = 0;
  virtual void correctionComponentsIs(const ReducedSeed * cc) = 0;

protected:
  UpdatedSeedAssembler() {}

  void setAssemblyPhase(PhaseRank ap) { assemblyPhase_ = ap; }
  void setUpdatedSeed(Seed * us) { updatedSeed_ = us; }
  void setPropagatedSeed(const Seed * ps) { propagatedSeed_ = ps; }
  void setCorrectionComponents(const ReducedSeed * cc) { correctionComponents_ = cc; }

private:
  PhaseRank assemblyPhase_;
  Seed::Ptr updatedSeed_;
  Seed::PtrConst propagatedSeed_;
  ReducedSeed::PtrConst correctionComponents_; 

  DISALLOW_COPY_AND_ASSIGN(UpdatedSeedAssembler);
};

} /* end namespace Pita */

#endif /* PITA_UPDATEDSEEDASSEMBLER_H */
