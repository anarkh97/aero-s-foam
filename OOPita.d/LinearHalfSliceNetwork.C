#include "LinearHalfSliceNetwork.h"
#include <Problems.d/DynamDescr.h>

#include "TrivialFullSlice.h"
#include "TrivialHalfSlice.h"

//#include "DynamStateReactorImpl.h"

namespace Pita {

LinearHalfSliceNetwork::LinearHalfSliceNetwork(CpuRank myCpu, HalfSliceMapping * m, SingleDomainDynamic<double> * pbDesc, HalfTimeSlice::Manager * halfSliceMgr) :
  HalfSliceNetwork(myCpu, m),
  probDesc_(pbDesc),
  halfSliceMgr_(halfSliceMgr)
{
  HalfSliceMapping::IteratorConst begin = mapping()->halfSlices(localCpuRank());
  for (HalfSliceMapping::IteratorConst it = begin; it; ++it) {
    // Do not instantiate inactive half-timeslices
    if (*it >= mapping()->firstInactiveHalfSlice())
      break;

    HalfSliceRank primalRank = *it;
    PrimalSliceRank primalSlice = primalRank.primalRank();
    
    HalfTimeSlice::Position pos;
    HalfSliceRank dualRank;
    if (primalSlice.headRank() == primalRank) {
      pos = HalfTimeSlice::HEAD;
      dualRank = primalRank + HalfSliceCount(1);
    } else {
      pos = HalfTimeSlice::TAIL;
      dualRank = primalRank - HalfSliceCount(1);
    }

    HalfTimeSlice::Ptr primalHs = halfSliceMgr_->instanceNew(HalfSliceId(primalRank, pos));
    DynamStateInput::Ptr seed = primalHs->seedInterface();
    DynamStateInterface::Ptr propagatedSeed = primalHs->propagatedSeedInterface();

    switch (pos) {
      case HalfTimeSlice::HEAD: {
        HalfSliceRank index = primalRank - HalfSliceCount(1);
        TrivialFullSliceTail::Ptr fullTail = new TrivialFullSliceTail(index);
        TrivialFullSliceHead::Ptr fullHead = new TrivialFullSliceHead(fullTail.ptr());
        fullTail_[index] = fullTail;
        fullHead_[index] = fullHead;

        //LocalDynamStateReactor::Ptr reactor = LocalDynamStateReactor::New(propagatedSeed.ptr(), fullTail->nextLeftPropagatedSeed().ptr());
      }
        break;
      case HalfTimeSlice::TAIL:
        break;
    }

    //HalfTimeSlice::Ptr dualHs = halfSliceMgr_->instanceNew(HalfSliceId(dualRank, pos));
  }

  for (FullTailMap::iterator it = fullTail_.begin(); it != fullTail_.end(); ++it) {
    // TODO Slots.... 
  }

}

} // end namespace Pita
