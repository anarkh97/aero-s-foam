#ifndef PITA_FINEINTEGRATORMANAGER_H
#define PITA_FINEINTEGRATORMANAGER_H

#include "Fwk.h"
#include "HalfTimeSlice.h"
#include "../DynamTimeIntegrator.h"

namespace Pita {

class FineIntegratorManager : public Fwk::PtrInterface<FineIntegratorManager> {
public:
  EXPORT_PTRINTERFACE_TYPES(FineIntegratorManager);
  
  typedef DynamTimeIntegrator IntegratorType;

  Seconds fineTimeStepSize() const { return fineTimeStepSize_; }
  IntegratorType * fineIntegrator(HalfTimeSlice::Direction direction) const;

protected:
  explicit FineIntegratorManager(Seconds fineTimeStep);

  virtual IntegratorType * createFineIntegrator(HalfTimeSlice::Direction direction) const = 0;

private:
  DISALLOW_COPY_AND_ASSIGN(FineIntegratorManager);

  Seconds fineTimeStepSize_;

  mutable IntegratorType::Ptr forwardFineIntegrator_;
  mutable IntegratorType::Ptr backwardFineIntegrator_;
};

} // end namespace Pita

#endif /* PITA_FINEINTEGRATORMANAGER_H */
