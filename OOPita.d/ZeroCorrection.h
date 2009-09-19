#ifndef PITA_ZEROCORRECTION_H
#define PITA_ZEROCORRECTION_H

#include "NamedTask.h"
#include "Seed.h"

namespace Pita {

class ZeroCorrection : public NamedTask {
public:
  EXPORT_PTRINTERFACE_TYPES(ZeroCorrection);

  const ReducedSeed * reducedCorrection() const { return reducedCorrection_.ptr(); }
  void reducedCorrectionIs(ReducedSeed * rc) { reducedCorrection_ = rc; }

  size_t reducedBasisSize() const { return reducedBasisSize_; }
  void reducedBasisSizeIs(size_t rbs) { reducedBasisSize_ = rbs; }

  // overriden
  virtual void iterationIs(IterationRank i) {
    reducedCorrection_->stateIs(Vector(reducedBasisSize(), 0.0));
    reducedCorrection_->statusIs(Seed::CONVERGED);
    reducedCorrection_->iterationIs(i);
  }

  static Ptr New() {
    return new ZeroCorrection();
  }

private:
  ReducedSeed::Ptr reducedCorrection_;
  size_t reducedBasisSize_;

  ZeroCorrection() :
    NamedTask("Zero Correction"),
    reducedBasisSize_(0)
  {}

  DISALLOW_COPY_AND_ASSIGN(ZeroCorrection);
};

} /* end namespace Pita */

#endif /* PITA_ZEROCORRECTION_H */
