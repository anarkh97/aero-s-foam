#ifndef PITA_LINEARDRIVER_H
#define PITA_LINEARDRIVER_H

#include "Fwk.h"
template <typename Scalar> class SingleDomainDynamic;

namespace Pita {

// Abstract base class to isolate main from the Pita module
class LinearDriver : public Fwk::PtrInterface<LinearDriver> {
public:
  EXPORT_PTRINTERFACE_TYPES(LinearDriver);

  SingleDomainDynamic<double> * probDesc() const { return probDesc_; }
  
  // Main routine
  virtual void solve() = 0;

protected:
  explicit LinearDriver(SingleDomainDynamic<double> * pbDesc) :
    probDesc_(pbDesc) {}

private:
  SingleDomainDynamic<double> * probDesc_;

  DISALLOW_COPY_AND_ASSIGN(LinearDriver);
};

} // end namespace Pita

extern Pita::LinearDriver::Ptr linearPitaDriverNew(SingleDomainDynamic<double> * pbDesc);
extern Pita::LinearDriver::Ptr linearReversiblePitaDriverNew(SingleDomainDynamic<double> * pbDesc);

#endif /* PITA_LINEARDRIVER_H */
