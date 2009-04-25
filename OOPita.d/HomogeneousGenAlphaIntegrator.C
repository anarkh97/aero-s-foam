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

} // end namespace Pita
