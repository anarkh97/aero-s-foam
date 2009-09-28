#ifndef PITA_HTS_NETWORKSTATICSTRATEGY_H
#define PITA_HTS_NETWORKSTATICSTRATEGY_H

#include "HalfSliceNetworkStrategy.h"

namespace Pita {

class HalfSliceNetworkStaticStrategy : public HalfSliceNetworkStrategy {
public:
  typedef Fwk::Ptr<HalfSliceNetworkStaticStrategy> Ptr;
  typedef Fwk::Ptr<const HalfSliceNetworkStaticStrategy> PtrConst;

  virtual HalfSliceRank workLoadShare(const Hts::SliceId & id) const;
  virtual SliceIteratorConst slice(HalfSliceRank workLoadId) const;

  static HalfSliceNetworkStaticStrategy::Ptr New() {
    return new HalfSliceNetworkStaticStrategy();
  }

protected:
  HalfSliceNetworkStaticStrategy() {}
};

} // end namespace Pita

#endif /* PITA_HTS_NETWORKSTATICSTRATEGY_H */
