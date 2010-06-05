#include "HomogeneousGenAlphaIntegrator.h"

namespace Pita {

HomogeneousGenAlphaIntegrator::HomogeneousGenAlphaIntegrator(
    LinearDynamOps::Manager * dOpsMgr,
    const GeneralizedAlphaParameter & param) :
  LinearGenAlphaIntegrator(dOpsMgr, param, HOMOGENEOUS)
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
  LinearGenAlphaIntegrator(dOpsMgr, param, NONHOMOGENEOUS)
{}

void
AffineGenAlphaIntegrator::externalForceStatusIs(AffineGenAlphaIntegrator::ExternalForceStatus efs) {
  zeroExternalForce();
  setExternalForceStatus(efs);
}

void
AffineGenAlphaIntegrator::computeExternalForce(
    Seconds forceEvalTime,
    SysState<VectorType> & currentState) { 
  if (externalForceStatus() == NONHOMOGENEOUS) {
    LinearGenAlphaIntegrator::computeExternalForce(forceEvalTime, currentState);
  }
}

} // end namespace Pita
