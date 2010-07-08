#ifndef PITA_HTS_PROPAGATIONDATASHARING_H
#define PITA_HTS_PROPAGATIONDATASHARING_H

#include "Fwk.h"
#include "Types.h"

#include "../NamedTask.h"

#include "SliceMapping.h"

#include "../DynamStateBasis.h"
#include "../DynamState.h"
#include "../Seed.h"

class Communicator;

#include "../SimpleBuffer.h"
#include <map>

#include "BasisCollectorImpl.h"

namespace Pita { namespace Hts {

class PropagationDataSharing : public NamedTask {
public:
  EXPORT_PTRINTERFACE_TYPES(PropagationDataSharing);

  // Execution
  virtual void iterationIs(IterationRank iter); // overriden

  // Parameters
  const SliceMapping * mapping() const { return mapping_.ptr(); }
  CpuRank localCpu() const { return localCpu_; }
  size_t vectorSize() const { return vectorSize_; }

  // Input
  BasisCollector * collector() { return collector_.ptr(); }

  // Output
  const DynamStateBasis * consolidatedBasis() const { return consolidatedBasis_.ptr(); }
   
  PropagationDataSharing(const SliceMapping * mapping, CpuRank localCpu, Communicator * timeComm, size_t vectorSize);

private:
  SliceMapping::PtrConst mapping_;
  CpuRank localCpu_;
  Communicator * timeComm_;
  size_t vectorSize_;

  BasisCollectorImpl::Ptr collector_;

  SimpleBuffer<double> buffer_;
  DynamStateBasis::Ptr consolidatedBasis_;
};

} /* end namespace Hts */ } /* end namespace Pita */

#endif /* PITA_HTS_PROPAGATIONDATASHARING_H */
