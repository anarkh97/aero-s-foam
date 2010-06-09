#ifndef PITA_LINEARDRIVER_H
#define PITA_LINEARDRIVER_H

#include "Fwk.h"
class SingleDomainDynamic;

namespace Pita {

class LinearDriver : public Fwk::PtrInterface<LinearDriver> {
public:
  EXPORT_PTRINTERFACE_TYPES(LinearDriver);

  // Main routine
  virtual void solve() = 0;

protected:
  LinearDriver() {}

private:
  DISALLOW_COPY_AND_ASSIGN(LinearDriver);
};

} // end namespace Pita

extern Pita::LinearDriver::Ptr linearPitaDriverNew(SingleDomainDynamic * pbDesc);

#endif /* PITA_LINEARDRIVER_H */
