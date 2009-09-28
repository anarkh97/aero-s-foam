#ifndef PITA_HTS_FULLTIMESLICE_HEADTAIL_H
#define PITA_HTS_FULLTIMESLICE_HEADTAIL_H

#include "Fwk.h"
#include "Types.h"
#include "../Seed.h"

namespace Pita { namespace Hts {

class FullTimeSliceHead : public Fwk::PtrInterface<FullTimeSliceHead> {
public:
  // Accessory types
  EXPORT_PTRINTERFACE_TYPES(FullTimeSliceHead);
  typedef Fwk::GenManagerInterface<FullTimeSliceHead *, HalfSliceRank> Manager;

  // Read accessors
  HalfSliceRank headHalfSlice() const { return headHalfSlice_; }
  HalfSliceRank tailHalfSlice() const { return tailHalfSlice_; }
  CpuRank tailCpu() const { return tailCpu_; }
  PhaseRank phase() const { return phase_; } 
  const Seed * updatedSeed() const { return updatedSeed_.ptr(); }
  const Seed * rightPropagatedSeed() const { return rightPropagatedSeed_.ptr(); }

  // Mutators
  virtual void tailCpuIs(CpuRank tc) { setTailCpu(tc); }
  virtual void phaseIs(PhaseRank p) { setPhase(p); }
  virtual void updatedSeedIs(const Seed * us) { setUpdatedSeed(us); }
  virtual void rightPropagatedSeedIs(const Seed * rps) { setRightPropagatedSeed(rps); }

protected:
  FullTimeSliceHead(HalfSliceRank headRank, CpuRank tailCpu);

  void setTailCpu(CpuRank tc) { tailCpu_ = tc; }
  void setPhase(PhaseRank p) { phase_ = p; }
  void setRightPropagatedSeed(const Seed * rps) { rightPropagatedSeed_ = rps; }
  void setUpdatedSeed(const Seed * us) { updatedSeed_ = us; }

private:
  HalfSliceRank headHalfSlice_;
  HalfSliceRank tailHalfSlice_;

  CpuRank tailCpu_;

  PhaseRank phase_;

  Fwk::Ptr<const Seed> updatedSeed_;
  Fwk::Ptr<const Seed> rightPropagatedSeed_;
 
  DISALLOW_COPY_AND_ASSIGN(FullTimeSliceHead); 
};


class FullTimeSliceTail : public Fwk::PtrInterface<FullTimeSliceTail> {
public:
  EXPORT_PTRINTERFACE_TYPES(FullTimeSliceTail);
  typedef Fwk::GenManagerInterface<FullTimeSliceTail *, HalfSliceRank> Manager;

  HalfSliceRank headHalfSlice() const { return headHalfSlice_; }
  HalfSliceRank tailHalfSlice() const { return tailHalfSlice_; }

  const Seed * nextUpdatedSeed() const { return nextUpdatedSeed_.ptr(); }
  Seed * nextUpdatedSeed() { return nextUpdatedSeed_.ptr(); }
  const Seed * nextLeftPropagatedSeed() const { return nextLeftPropagatedSeed_.ptr(); }

  CpuRank headCpu() const { return headCpu_; }

  PhaseRank phase() const { return phase_; }

  virtual void nextUpdatedSeedIs(Seed * seed) { setNextUpdatedSeed(seed); }
  virtual void nextLeftPropagatedSeedIs(const Seed * seed) { setNextLeftPropagatedSeed(seed); }
  
  virtual void headCpuIs(CpuRank hc) { setHeadCpu(hc); }
  
  virtual void phaseIs(PhaseRank p) { setPhase(p); }

protected:
  FullTimeSliceTail(HalfSliceRank tailRank, CpuRank headCpu);
  
  void setNextLeftPropagatedSeed(const Seed * nlps) { nextLeftPropagatedSeed_ = nlps; }
  void setNextUpdatedSeed(Seed * nus) { nextUpdatedSeed_ = nus; }
  void setHeadCpu(CpuRank hc) { headCpu_ = hc; }
  void setPhase(PhaseRank p) { phase_ = p; }

private:
  HalfSliceRank headHalfSlice_;
  HalfSliceRank tailHalfSlice_;
  
  CpuRank headCpu_;
 
  Fwk::Ptr<const Seed> nextLeftPropagatedSeed_;
  Fwk::Ptr<Seed> nextUpdatedSeed_;
 
  PhaseRank phase_;

  DISALLOW_COPY_AND_ASSIGN(FullTimeSliceTail);
};

} /* end namespace Hts */ } /* end namespace Pita */

#endif /* PITA_HTS_FULLTIMESLICE_HEADTAIL_H */
