#ifndef PITA_INTEGRATORPROPAGATOR_H
#define PITA_INTEGRATORPROPAGATOR_H

#include "Fwk.h"
#include "Types.h"
#include "DynamPropagator.h" 
#include "DynamTimeIntegrator.h"

namespace Pita {

class IntegratorPropagator : public DynamPropagator {
public:
  EXPORT_PTRINTERFACE_TYPES(IntegratorPropagator);

  DynamTimeIntegrator * integrator() const { return integrator_.ptr(); }
  Seconds initialTime() const { return initialTime_; }
  TimeStepCount timeStepCount() const { return timeStepCount_; }
  Seconds timeStepSize() const { return timeStepSize_; }
  
  virtual void initialStateIs(const DynamState & initialState);
  void integratorIs(DynamTimeIntegrator * integrator) { integrator_ = integrator; }
  void initialTimeIs(Seconds t0) { initialTime_ = t0; }
  void timeStepCountIs(TimeStepCount timeStepCount) { timeStepCount_ = timeStepCount; }
  void timeStepSizeIs(Seconds tss) { timeStepSize_ = tss; }

  // HACK
  SliceRank sliceRank() const { return sliceRank_; }
  void sliceRankIs(SliceRank slice) { sliceRank_ = slice; }
  
  static IntegratorPropagator::Ptr New(DynamTimeIntegrator * integrator) {
    return new IntegratorPropagator(integrator);
  }
  
protected:
  explicit IntegratorPropagator(DynamTimeIntegrator * integrator);

private:
  Fwk::Ptr<DynamTimeIntegrator> integrator_;
  Seconds initialTime_;
  TimeStepCount timeStepCount_;
  Seconds timeStepSize_;

  // HACK
  SliceRank sliceRank_;
};
  
} // end namespace Pita

#endif /* PITA_INTEGRATORPROPAGATOR_H */
