#ifndef PITA_HOMOGENEOUSGENALPHAINTEGRATOR_H
#define PITA_HOMOGENEOUSGENALPHAINTEGRATOR_H

#include "LinearGenAlphaIntegrator.h"

namespace Pita {

class HomogeneousGenAlphaIntegrator : public LinearGenAlphaIntegrator {
public:
  EXPORT_PTRINTERFACE_TYPES(HomogeneousGenAlphaIntegrator);

  HomogeneousGenAlphaIntegrator(LinearDynamOps::Manager * dOpsMgr, const GeneralizedAlphaParameter & param);

protected:
  // Overriden
  virtual void computeExternalForce(Seconds forceEvalTime, SysState<VectorType> & currentState);
};

} // end namespace Pita

#endif /* PITA_HOMOGENEOUSGENALPHAINTEGRATOR_H */
