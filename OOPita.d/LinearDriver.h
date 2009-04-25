#ifndef PITA_LINEARDRIVER_H
#define PITA_LINEARDRIVER_H

#include "Fwk.h"
template <typename Scalar> class SingleDomainDynamic;

namespace Pita {

class LinearDriver : public Fwk::PtrInterface<LinearDriver> {
public:
  typedef Fwk::Ptr<LinearDriver> Ptr;
  typedef Fwk::Ptr<const LinearDriver> PtrConst;

  // Main routine
  virtual void solve() = 0;
};

} // end namespace Pita

extern Pita::LinearDriver::Ptr linearPitaDriverNew(SingleDomainDynamic<double> * pbDesc);

#endif /* PITA_LINEARDRIVER_H */
