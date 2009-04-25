#ifndef PITA_LINEARIZEDPROPAGATOR_H
#define PITA_LINEARIZEDPROPAGATOR_H

#include "DynamPropagator.h"

namespace Pita {

class PitaNonLinDynamic;
class NlDynamTimeIntegrator;

class LinearizedPropagator : public DynamPropagator {
public:
  typedef Fwk::Ptr<LinearizedPropagator> Ptr;
  typedef Fwk::Ptr<const LinearizedPropagator> PtrConst;
  
  virtual void initialStateIs(const DynamState & initialState);

  static LinearizedPropagator::Ptr New(const NlDynamTimeIntegrator * integrator) {
    return new LinearizedPropagator(integrator);
  }

protected:
  explicit LinearizedPropagator(const NlDynamTimeIntegrator * integrator);

private: 
  Fwk::Ptr<const NlDynamTimeIntegrator> integrator_;
  const PitaNonLinDynamic * probDesc_;
  GenVector<double> temp_;
  GenVector<double> midTimeDisp_;
};  
  
} // end namespace Pita

#endif /* PITA_LINEARIZEDPROPAGATOR_H */
