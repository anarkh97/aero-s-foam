#ifndef PITA_LINEARHALFSLICENETWORK_H
#define PITA_LINEARHALFSLICENETWORK_H

#include "HalfSliceNetwork.h"
#include "HalfSliceSlot.h"
#include "FullSliceSlot.h"
#include <map>

template <typename Scalar> class SingleDomainDynamic;

namespace Pita {

class LinearHalfSliceNetwork : public HalfSliceNetwork {
public:
  typedef Fwk::Ptr<LinearHalfSliceNetwork> Ptr;
  typedef Fwk::Ptr<const LinearHalfSliceNetwork> PtrConst;

  const HalfTimeSlice::Manager * halfSliceManager() const { return halfSliceMgr_.ptr(); }
  HalfTimeSlice::Manager * halfSliceManager() { return halfSliceMgr_.ptr(); }

  static LinearHalfSliceNetwork::Ptr New(CpuRank myCpu,
                                         HalfSliceMapping * sMap,
                                         SingleDomainDynamic<double> * pbDesc,
                                         HalfTimeSlice::Manager * halfSliceMgr) {
    return new LinearHalfSliceNetwork(myCpu, sMap, pbDesc, halfSliceMgr);
  }

protected:
  LinearHalfSliceNetwork(CpuRank myCpu,
                         HalfSliceMapping * m,
                         SingleDomainDynamic<double> * pbDesc,
                         HalfTimeSlice::Manager * halfSliceMgr);

private:
  SingleDomainDynamic<double> * probDesc_;
  HalfTimeSlice::Manager::Ptr halfSliceMgr_;

  typedef std::map<HalfSliceId, HalfSliceSlot> HalfSliceSlotMap;
  HalfSliceSlotMap halfSliceSlot_;

  typedef std::map<HalfSliceRank, FullTimeSliceHead::Ptr> FullHeadMap;
  FullHeadMap fullHead_;

  typedef std::map<HalfSliceRank, FullTimeSliceTail::Ptr> FullTailMap;
  FullTailMap fullTail_;
};

} // end namespace Pita

#endif /* PITA_LINEARHALFSLICENETWORK_H */
