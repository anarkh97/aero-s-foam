#ifndef PITA_HTS_LINEARDRIVERIMPL_H
#define PITA_HTS_LINEARDRIVERIMPL_H

#include "../LinearDriver.h"

Pita::LinearDriver::Ptr linearPitaDriverNew(SingleDomainDynamic<double> * pbDesc);

namespace Pita { namespace Hts {

class LinearDriverImpl : public LinearDriver {
public:
  typedef Fwk::Ptr<LinearDriverImpl> Ptr;
  typedef Fwk::Ptr<const LinearDriverImpl> PtrConst;

  virtual void solve();

  SingleDomainDynamic<double> * probDesc() const { return probDesc_; }

  static LinearDriverImpl::Ptr New(SingleDomainDynamic<double> * pbDesc) {
    return new LinearDriverImpl(pbDesc);
  }

protected:
  explicit LinearDriverImpl(SingleDomainDynamic<double> * pbDesc);

private:
  SingleDomainDynamic<double> * probDesc_;
};

} /* end namespace namespace Hts */ } /* end namespace Pita */

#endif /* PITA_HTS_LINEARDRIVERIMPL_H */
