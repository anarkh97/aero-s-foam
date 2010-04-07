#ifndef PITA_HTS_CORRECTIONTIMESLICE_H
#define PITA_HTS_CORRECTIONTIMESLICE_H

#include "Fwk.h"
#include "Types.h"
#include "../NamedTask.h"
#include "../Seed.h"

namespace Pita { namespace Hts {

class CorrectionTimeSlice : public NamedTask {
public:
  EXPORT_PTRINTERFACE_TYPES(CorrectionTimeSlice);
  typedef Fwk::GenManagerInterface<CorrectionTimeSlice *, HalfSliceRank> Manager;
  
  // Id
  HalfSliceRank headHalfSlice() const { return headHalfSlice_; }
  HalfSliceRank tailHalfSlice() const { return tailHalfSlice_; }

  // Inputs  
  const Seed * predictedSeed() const { return predictedSeed_.ptr(); }
  const Seed * actualSeed() const { return actualSeed_.ptr(); } 
  const Seed * correction() const { return correction_.ptr(); }
 
  virtual void predictedSeedIs(const Seed * ps) { setPredictedSeed(ps); }
  virtual void actualSeedIs(const Seed * as) { setActualSeed(as); }
  virtual void correctionIs(const Seed * c) { setCorrection(c); }

  // Outputs
  Seed * jump() const { return jump_.ptr(); }
  Seed * nextCorrection() const { return nextCorrection_.ptr(); }
 
  virtual void jumpIs(Seed * j) { setJump(j); } 
  virtual void nextCorrectionIs(Seed * nc) { setNextCorrection(nc); }
  
protected:
  explicit CorrectionTimeSlice(HalfSliceRank headRank) :
    NamedTask(String("CorrectionTimeSlice") + toString(headRank)),
    headHalfSlice_(headRank),
    tailHalfSlice_(headRank + HalfSliceCount(1))
  {}
  
  void setPredictedSeed(const Seed * ps) { predictedSeed_ = ps; }
  void setActualSeed(const Seed * as) { actualSeed_ = as; }
  void setCorrection(const Seed * c) { correction_ = c; }
  void setJump(Seed * j) { jump_ = j; }
  void setNextCorrection(Seed * nc) { nextCorrection_ = nc; }
  
private:
  HalfSliceRank headHalfSlice_;
  HalfSliceRank tailHalfSlice_;
  
  Seed::PtrConst predictedSeed_;
  Seed::PtrConst actualSeed_;
  Seed::PtrConst correction_;
  Seed::Ptr jump_;
  Seed::Ptr nextCorrection_;

  DISALLOW_COPY_AND_ASSIGN(CorrectionTimeSlice);
};
  
} /* end namespace Hts */ } /* end namespace Pita */



#endif /* PITA_HTS_CORRECTIONTIMESLICE_H */
