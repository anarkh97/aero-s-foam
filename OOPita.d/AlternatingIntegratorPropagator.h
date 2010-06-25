#ifndef PITA_ALTERNATINGINTEGRATORPROPAGATOR_H
#define PITA_ALTERNATINGINTEGRATORPROPAGATOR_H

#include "Fwk.h"
#include "IntegratorPropagator.h"
#include "ProjectorPropagator.h"
#include "DynamStatePlainBasis.h"

namespace Pita {

class AlternatingIntegratorPropagator : public IntegratorPropagator {
public:
  EXPORT_PTRINTERFACE_TYPES(AlternatingIntegratorPropagator);

  virtual void initialStateIs(const DynamState & initialState);

  // New members
  const ProjectorPropagator * secondaryPropagator() const { return secondaryPropagator_.ptr(); }
  void secondaryPropagatorIs(ProjectorPropagator::Ptr secondary) { secondaryPropagator_ = secondary; }
  
  static AlternatingIntegratorPropagator::Ptr New(DynamTimeIntegrator * integrator) {
    return new AlternatingIntegratorPropagator(integrator);
  }
  
protected:
  explicit AlternatingIntegratorPropagator(DynamTimeIntegrator * integrator);
  
private:
  ProjectorPropagator::Ptr secondaryPropagator_;
};

} // end namespace Pita

#endif /* PITA_ALTERNATINGINTEGRATORPROPAGATOR_H */
