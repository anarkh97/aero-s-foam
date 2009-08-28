#ifndef PITA_REDUCEDFULLTIMESLICE_H
#define PITA_REDUCEDFULLTIMESLICE_H

#include "Fwk.h"
#include "Types.h"
#include "Seed.h"

namespace Pita {

class ReducedFullTimeSlice : public Fwk::PtrInterface<ReducedFullTimeSlice> {
public:
  // Accessory types
  EXPORT_PTRINTERFACE_TYPES(ReducedFullTimeSlice);
  typedef Fwk::GenManagerInterface<ReducedFullTimeSlice *, HalfSliceRank> Manager;

  // Id
  HalfSliceRank headHalfSlice() const { return headHalfSlice_; }
  HalfSliceRank tailHalfSlice() const { return tailHalfSlice_; }
  
  // Seeds
  const ReducedSeed * jump() const { return jump_.ptr(); }
  const ReducedSeed * correction() const { return correction_.ptr(); }
  const ReducedSeed * nextCorrection() const { return nextCorrection_.ptr(); }
  ReducedSeed * nextCorrection() { return nextCorrection_.ptr(); }

  virtual void jumpIs(const ReducedSeed * j) { setJump(j); }
  virtual void correctionIs(const ReducedSeed * c) { setCorrection(c); }
  virtual void nextCorrectionIs(ReducedSeed * nc) { setNextCorrection(nc); }

  // Iteration control
  IterationRank iteration() const { return iteration_; }
  virtual void iterationIs(IterationRank i) = 0;

protected:
  explicit ReducedFullTimeSlice(HalfSliceRank headRank) :
    headHalfSlice_(headRank),
    tailHalfSlice_(headRank + HalfSliceCount(1))
  {}

  void setJump(const ReducedSeed * j) { jump_ = j; }
  void setCorrection(const ReducedSeed * c) { correction_ = c; }
  void setNextCorrection(ReducedSeed * c) { nextCorrection_ = c; }
  void setIteration(IterationRank i) { iteration_ = i; }

private:
  HalfSliceRank headHalfSlice_;
  HalfSliceRank tailHalfSlice_;

  Fwk::Ptr<const ReducedSeed> jump_;
  Fwk::Ptr<const ReducedSeed> correction_;
  Fwk::Ptr<ReducedSeed> nextCorrection_;

  IterationRank iteration_;

  DISALLOW_COPY_AND_ASSIGN(ReducedFullTimeSlice); 
};

} // end namespace Pita

#endif /* PITA_REDUCEDFULLTIMESLICE_H */
