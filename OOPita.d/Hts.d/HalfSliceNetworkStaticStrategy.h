#ifndef PITA_HALFSLICENETWORKSTATICSTRATEGY_H
#define PITA_HALFSLICENETWORKSTATICSTRATEGY_H

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

#endif /* PITA_HALFSLICENETWORKSTATICSTRATEGY_H */
