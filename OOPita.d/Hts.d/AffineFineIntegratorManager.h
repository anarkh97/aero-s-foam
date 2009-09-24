#ifndef PITA_HTS_AFFINEFINEINTEGRATORMANAGER_H
#define PITA_HTS_AFFINEFINEINTEGRATORMANAGER_H

#include "LinearFineIntegratorManager.h"

#include "../HomogeneousGenAlphaIntegrator.h"

#include "../Activity.h"
#include "HalfSliceSchedule.h"

#include <list>

namespace Pita { namespace Hts {

class AffineFineIntegratorManager : public LinearFineIntegratorManager<AffineGenAlphaIntegrator> {
public:
  EXPORT_PTRINTERFACE_TYPES(AffineFineIntegratorManager);

  const Schedule * schedule() const { return schedule_.ptr(); }

  static Ptr New(LinearDynamOps::Manager * dom, const GeneralizedAlphaParameter & fp, const Schedule * schedule) {
    return new AffineFineIntegratorManager(dom, fp, schedule);
  }
  
protected:
  AffineFineIntegratorManager(LinearDynamOps::Manager * dom, const GeneralizedAlphaParameter & fp, const Schedule * schedule);
  
  class AffineForceReactor : public Activity::Notifiee {
  public:
    EXPORT_PTRINTERFACE_TYPES(AffineForceReactor);

    virtual void onStatus();

    static Ptr New(Activity * notifier, AffineGenAlphaIntegrator * target) {
      return new AffineForceReactor(notifier, target);
    }
    
  protected:
    AffineForceReactor(Activity * notifier, AffineGenAlphaIntegrator * target); 

  private:
    AffineGenAlphaIntegrator::Ptr target_;   
  };

  virtual AffineGenAlphaIntegrator * createFineIntegrator(HalfTimeSlice::Direction direction) const; // Overriden

private:
  Schedule::PtrConst schedule_;

  mutable std::list<AffineForceReactor::Ptr> reactor_;
};

} /* end namespace Hts */ } /* end namespace Pita */

#endif /* PITA_HTS_AFFINEFINEINTEGRATORMANAGER_H */
