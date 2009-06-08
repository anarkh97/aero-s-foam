#include "HomogeneousGenAlphaIntegrator.h"

namespace Pita {

HomogeneousGenAlphaIntegrator::HomogeneousGenAlphaIntegrator(
    LinearDynamOps::Manager * dOpsMgr,
    const GeneralizedAlphaParameter & param) :
  LinearGenAlphaIntegrator(dOpsMgr, param)
{}

void
HomogeneousGenAlphaIntegrator::computeExternalForce(
    Seconds forceEvalTime,
    SysState<VectorType> & currentState) { 
  // Do nothing
}

AffineGenAlphaIntegrator::AffineGenAlphaIntegrator(
    LinearDynamOps::Manager * dOpsMgr,
    const GeneralizedAlphaParameter & param) :
  LinearGenAlphaIntegrator(dOpsMgr, param),
  externalForceFlag_(true)
{}

void
AffineGenAlphaIntegrator::externalForceFlagIs(bool eff) {
  externalForceFlag_ = eff;
  zeroExternalForce();
}

void
AffineGenAlphaIntegrator::computeExternalForce(
    Seconds forceEvalTime,
    SysState<VectorType> & currentState) { 
  if (externalForceFlag()) {
    LinearGenAlphaIntegrator::computeExternalForce(forceEvalTime, currentState);
  }
}

} // end namespace Pita
