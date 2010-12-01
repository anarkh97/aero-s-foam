#ifndef PITA_HTS_GLOBALSTATESHARING_H
#define PITA_HTS_GLOBALSTATESHARING_H

#include "Fwk.h"
#include "Types.h"

#include "../NamedTask.h"

#include "SliceMapping.h"

#include "../DynamStateBasis.h"
#include "../DynamState.h"
#include "../Seed.h"

#include "LocalNetworkImpl.h"

class Communicator;

#include "../SimpleBuffer.h"
#include <map>
#include <queue>

namespace Pita { namespace Hts {

class GlobalStateSharing : public NamedTask {
public:
  EXPORT_PTRINTERFACE_TYPES(GlobalStateSharing);

  // Parameters
  CpuRank localCpu() const;
  size_t vectorSize() const { return vectorSize_; }

  // Input
  typedef LocalNetworkImpl::SeedGetter<DynamState> SeedGetter;
  const SeedGetter & seedGetter() const { return seedGetter_; }
  void seedGetterIs(SeedGetter sg) { seedGetter_ = sg; }

  // Output
  const DynamStateBasis * consolidatedBasis() const { return consolidatedBasis_.ptr(); }
  
  // Setup exchange parameters
  virtual void mappingIs(const SliceMapping & m);
  
  // Execution
  virtual void iterationIs(IterationRank iter); // overriden

  GlobalStateSharing(Communicator * timeComm, size_t vectorSize);

private:
  Communicator * timeComm_;
  size_t vectorSize_;

  SeedGetter seedGetter_;

  std::queue<Seed::PtrConst> localStates_;
  int stateCount_;
  SimpleBuffer<double> buffer_;
  SimpleBuffer<int> bufferCounts_;
  SimpleBuffer<int> bufferStrides_; 

  DynamStateBasis::Ptr consolidatedBasis_;
};

} /* end namespace Hts */ } /* end namespace Pita */

#endif /* PITA_HTS_GLOBALSTATESHARING_H */
