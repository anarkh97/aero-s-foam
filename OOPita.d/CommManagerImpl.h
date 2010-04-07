#ifndef PITA_COMMMANAGERIMPL_H
#define PITA_COMMMANAGERIMPL_H

#include "CommManager.h"
#include "SimpleBuffer.h"
#include <algorithm>

class Communicator;

namespace Pita {

class CommManagerImpl : public CommManager {
public:
  typedef Fwk::Ptr<CommManagerImpl> Ptr;
  typedef Fwk::Ptr<const CommManagerImpl> PtrConst;

  virtual void localBroadcastBasisIs(DynamStateBasis::Ptr basis); 
  
  virtual void correctionVectorSizeIs(size_t s);
  virtual void nextCorrectionIs(const DynamState & state);
  
  static CommManagerImpl::Ptr New(TimeSliceMapping * sliceMap, Communicator * timeComm);
  
protected:
  CommManagerImpl(TimeSliceMapping * sliceMap, CpuRank myCpuRank, CpuRank nextCpuRank, CpuRank previousCpuRank, Communicator * timeComm);

private:
  Communicator * timeCommunicator_;
  SimpleBuffer<double> gBuffer_;
  SimpleBuffer<double> sBuffer_;
  SimpleBuffer<double> rBuffer_;
  SimpleBuffer<int> mpiParameters_;

  void receiveCorrection();
  void sendCorrection();
};

} // end namespace Pita

#endif /* PITA_COMMMANAGERIMPL_H */
