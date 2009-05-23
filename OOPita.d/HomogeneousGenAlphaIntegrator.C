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

HackedGenAlphaIntegrator::HackedGenAlphaIntegrator(
    LinearDynamOps::Manager * dOpsMgr,
    const GeneralizedAlphaParameter & param) :
  LinearGenAlphaIntegrator(dOpsMgr, param),
  externalForceFlag(true)
{}

void
HackedGenAlphaIntegrator::computeExternalForce(
    Seconds forceEvalTime,
    SysState<VectorType> & currentState) { 
  if (externalForceFlag) {
    LinearGenAlphaIntegrator::computeExternalForce(forceEvalTime, currentState);
    externalForceFlag = false;
  }
}

} // end namespace Pita
