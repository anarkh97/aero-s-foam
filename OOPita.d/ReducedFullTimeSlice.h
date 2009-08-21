#ifndef PITA_REDUCEDFULLTIMESLICE_H
#define PITA_REDUCEDFULLTIMESLICE_H

#include "Fwk.h"
#include "Types.h"
#include "Seed.h"

namespace Pita {

class ReducedFullTimeSliceHead : public Fwk::PtrInterface<ReducedFullTimeSliceHead> {
public:
  // Accessory types
  EXPORT_PTRINTERFACE_TYPES(ReducedFullTimeSliceHead);
  typedef Fwk::GenManagerInterface<ReducedFullTimeSliceHead *, HalfSliceRank> Manager;

  // Read accessors
  HalfSliceRank headHalfSlice() const { return headHalfSlice_; }
  HalfSliceRank tailHalfSlice() const { return tailHalfSlice_; }
  CpuRank tailCpu() const { return tailCpu_; }
  PhaseRank phase() const { return phase_; } 
  const Seed * rightPropagatedSeed() const { return rightPropagatedSeed_.ptr(); }
  const Seed * leftPropagatedSeed() const { return leftPropagatedSeed_.ptr(); }
  const Seed * updatedSeed() const { return updatedSeed_.ptr(); }
  Seed * updatedSeed() { return updatedSeed_.ptr(); }
  const ReducedSeed * correction() const { return correction_.ptr(); }

  // Mutators
  virtual void tailCpuIs(CpuRank tc) { setTailCpu(tc); }
  virtual void phaseIs(PhaseRank p) { setPhase(p); }
  virtual void rightPropagatedSeedIs(const Seed * rps) { setRightPropagatedSeed(rps); }
  virtual void leftPropagatedSeedIs(const Seed * lps) { setLeftPropagatedSeed(lps); }
  virtual void updatedSeedIs(Seed * us) { setUpdatedSeed(us); }
  virtual void correctionIs(const ReducedSeed * c) { setCorrection(c); }

protected:
  ReducedFullTimeSliceHead(HalfSliceRank headRank, CpuRank tailCpu);

  void setTailCpu(CpuRank tc) { tailCpu_ = tc; }
  void setPhase(PhaseRank p) { phase_ = p; }
  void setRightPropagatedSeed(const Seed * rps) { rightPropagatedSeed_ = rps; }
  void setLeftPropagatedSeed(const Seed * lps) { leftPropagatedSeed_ = lps; }
  void setUpdatedSeed(Seed * us) { updatedSeed_ = us; }
  void setCorrection(const ReducedSeed * c) { correction_ = c; }

private:
  HalfSliceRank headHalfSlice_;
  HalfSliceRank tailHalfSlice_;

  CpuRank tailCpu_;

  PhaseRank phase_;

  Fwk::Ptr<const Seed> rightPropagatedSeed_;
  Fwk::Ptr<const Seed> leftPropagatedSeed_;
  Fwk::Ptr<Seed> updatedSeed_;
  Fwk::Ptr<const ReducedSeed> correction_;
 
  DISALLOW_COPY_AND_ASSIGN(ReducedFullTimeSliceHead); 
};


class ReducedFullTimeSliceTail : public Fwk::PtrInterface<ReducedFullTimeSliceTail> {
public:
  EXPORT_PTRINTERFACE_TYPES(ReducedFullTimeSliceTail);
  typedef Fwk::GenManagerInterface<ReducedFullTimeSliceTail *, HalfSliceRank> Manager;

  HalfSliceRank headHalfSlice() const { return headHalfSlice_; }
  HalfSliceRank tailHalfSlice() const { return tailHalfSlice_; }

  const Seed * nextLeftPropagatedSeed() const { return nextLeftPropagatedSeed_.ptr(); }
  const Seed * nextUpdatedSeed() const { return nextUpdatedSeed_.ptr(); }
  Seed * nextUpdatedSeed() { return nextUpdatedSeed_.ptr(); }
  const ReducedSeed * nextCorrection() const { return nextCorrection_.ptr(); }
  ReducedSeed * nextCorrection() { return nextCorrection_.ptr(); }

  CpuRank headCpu() const { return headCpu_; }

  PhaseRank phase() const { return phase_; }

  virtual void nextLeftPropagatedSeedIs(const Seed * nlps) { setNextLeftPropagatedSeed(nlps); }
  virtual void nextUpdatedSeedIs(Seed * nus) { setNextUpdatedSeed(nus); }
  virtual void nextCorrectionIs(ReducedSeed * nc) { setNextCorrection(nc); }

  virtual void headCpuIs(CpuRank hc) { setHeadCpu(hc); }
  
  virtual void phaseIs(PhaseRank p) { setPhase(p); }

protected:
  ReducedFullTimeSliceTail(HalfSliceRank tailRank, CpuRank headCpu);
  
  void setNextLeftPropagatedSeed(const Seed * nlps) { nextLeftPropagatedSeed_ = nlps; }
  void setNextUpdatedSeed(Seed * nus) { nextUpdatedSeed_ = nus; }
  void setNextCorrection(ReducedSeed * nc) { nextCorrection_ = nc; }
  void setHeadCpu(CpuRank hc) { headCpu_ = hc; }
  void setPhase(PhaseRank p) { phase_ = p; }

private:
  HalfSliceRank headHalfSlice_;
  HalfSliceRank tailHalfSlice_;
  
  CpuRank headCpu_;
 
  Fwk::Ptr<const Seed> nextLeftPropagatedSeed_;
  Fwk::Ptr<Seed> nextUpdatedSeed_;
  Fwk::Ptr<ReducedSeed> nextCorrection_;
 
  PhaseRank phase_;

  DISALLOW_COPY_AND_ASSIGN(ReducedFullTimeSliceTail);
};

} // end namespace Pita

#endif /* PITA_REDUCEDFULLTIMESLICE_H */
