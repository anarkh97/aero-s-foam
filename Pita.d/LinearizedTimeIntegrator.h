#ifndef _LINEARIZED_TIME_INTEGRATOR_H_
#define _LINEARIZED_TIME_INTEGRATOR_H_

class PitaNonLinDynamic;
template <typename Scalar> class DynamStateSet;

#include <Math.d/Vector.h>
#include <Pita.d/DynamState.h>

class LinearizedTimeIntegrator
{
public:
  typedef double Scalar;
  typedef GenVector<Scalar> VecType;
  explicit LinearizedTimeIntegrator(PitaNonLinDynamic &);
  void stepLinearizedIntegrate(DynamStateSet<double> &) const;

private:
  PitaNonLinDynamic & probDesc;
  double delta;
  double oneOverDelta;
  mutable VecType midTimeDisp;
  mutable VecType temp;
  mutable DynamState<Scalar> newState;
};

#endif
