#ifndef PITA_LINEARGENALPHAINTEGRATOR_H
#define PITA_LINEARGENALPHAINTEGRATOR_H

#include "AffineDynamTimeIntegrator.h"
#include "LinearDynamOps.h"
template <typename Scalar> class SingleDomainDynamic;
template <typename VectorType> class SysState;

namespace Pita {

class LinearGenAlphaIntegrator : public AffineDynamTimeIntegrator {
public:
  EXPORT_PTRINTERFACE_TYPES(LinearGenAlphaIntegrator);

  // Type aliases for (future) templating
  typedef double Scalar;
  typedef ::GenVector<Scalar> VectorType;
  typedef SingleDomainDynamic<Scalar> ProblemDescriptor;

  // Overriden methods
  virtual void timeStepSizeIs(Seconds dt);
  virtual void initialConditionIs(const DynamState & initialState, Seconds initialTime = Seconds(0.0));
  virtual void currentTimeInc(Seconds timeIncrement);
  virtual void timeStepCountInc(TimeStepCount steps = TimeStepCount(1));
  virtual void externalForceStatusIs(ExternalForceStatus efs);
  
  // Additional methods
  double rhoInfinity() const { return rhoInfinity_; }
  void rhoInfinityIs(double r);

  const LinearDynamOps * dynamOps() const { return dynamOps_.ptr(); }
  const VectorType & currentAcceleration() const { return currentAcceleration_; }
  void currentAccelerationIs(const VectorType & accel) { currentAcceleration_ = accel; }
  const VectorType & externalForce() const { return externalForce_; }
  const VectorType & previousVelocity() const { return previousVelocity_; }
  const VectorType & aeroForce() const { return aeroForce_; } 

  static LinearGenAlphaIntegrator::Ptr New(LinearDynamOps::Manager * dOpsMgr, const GeneralizedAlphaParameter & param) {
    return new LinearGenAlphaIntegrator(dOpsMgr, param, NONHOMOGENEOUS);
  }

protected:
  LinearGenAlphaIntegrator(LinearDynamOps::Manager * dOpsMgr, const GeneralizedAlphaParameter & param, ExternalForceStatus efs);
  
  const ProblemDescriptor * probDesc() const { return probDesc_; }
  ProblemDescriptor * probDesc() { return probDesc_; }

  void zeroExternalForce() { externalForce_.zero(); }
    
  virtual void computeExternalForce(Seconds forceEvalTime, SysState<VectorType> & currentState); 

private:
  void getDynamOps(const GeneralizedAlphaParameter & param);
  void integrate(TimeStepCount s);

  LinearDynamOps::Manager::Ptr dynamOpsManager_;
  ProblemDescriptor * probDesc_;
  LinearDynamOps::Ptr dynamOps_;

  double rhoInfinity_;
  double beta_, gamma_, alphaf_, alpham_;

  VectorType constForce_;
  VectorType externalForce_;
  VectorType aeroForce_;
  VectorType currentDisplacement_;
  VectorType currentVelocity_;
  VectorType currentAcceleration_;
  VectorType nextDisplacement_;
  VectorType nextVelocity_;
  VectorType nextAcceleration_;
  VectorType previousVelocity_;
  VectorType rhs_;
  VectorType temp_;
  VectorType temp2_;
};


} // end namespace Pita

#endif /* PITA_LINEARGENALPHAINTEGRATOR_H */
