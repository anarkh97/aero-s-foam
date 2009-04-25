#ifndef PITA_LINEARFINEINTEGRATORMANAGER_H
#define PITA_LINEARFINEINTEGRATORMANAGER_H

#include "FineIntegratorManager.h"

#include "LinearDynamOps.h"
#include "LinearGenAlphaIntegrator.h"

namespace Pita {

class LinearFineIntegratorManager : public FineIntegratorManager {
public:
  EXPORT_PTRINTERFACE_TYPES(LinearFineIntegratorManager);

  LinearDynamOps::Manager * dynamOpsManager() const { return dynamOpsManager_.ptr(); }
  const GeneralizedAlphaParameter & parameter(HalfTimeSlice::Direction direction) const;

  static Ptr New(LinearDynamOps::Manager * dynOpsMgr, const GeneralizedAlphaParameter & forwardParam) {
    return new LinearFineIntegratorManager(dynOpsMgr, forwardParam);
  }

protected:
  LinearFineIntegratorManager(LinearDynamOps::Manager * dom, const GeneralizedAlphaParameter & fp);

  virtual LinearGenAlphaIntegrator * createFineIntegrator(HalfTimeSlice::Direction direction) const; // Overriden

private:
  LinearDynamOps::Manager::Ptr dynamOpsManager_;
  GeneralizedAlphaParameter forwardParameter_;
  GeneralizedAlphaParameter backwardParameter_;
};

} // end namespace Pita

#endif /* PITA_LINEARFINEINTEGRATORMANAGER_H */
