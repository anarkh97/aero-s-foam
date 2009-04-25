#ifndef PITA_NLDRIVERIMPL_H
#define PITA_NLDRIVERIMPL_H

#include "NlDriver.h"
#include "PitaNonLinDynam.h"
#include "CommManager.h"
#include "TimeSliceMapping.h"
#include "TimeSliceNetwork.h"
#include "SeedInitializer.h"
#include "NlDynamPostProcessor.h"

Pita::NlDriver::Ptr nlPitaDriverNew(Pita::PitaNonLinDynamic * problemDescriptor);

namespace Pita {

class NlDriverImpl : public NlDriver {
public:
  typedef Fwk::Ptr<NlDriverImpl> Ptr;
  typedef Fwk::Ptr<const NlDriverImpl> PtrConst;

  virtual void solve();

  explicit NlDriverImpl(PitaNonLinDynamic * problemDescriptor);

protected:
  size_t vectorSize() const { return probDesc()->solVecInfo(); }
  CpuRank myCpuRank() const { return CpuRank(timeCommunicator_->myID()); }
  
private:
  Communicator * timeCommunicator_;

  void duplicateFiles(const TimeSliceMapping *);
};

} // end namespace Pita

#endif /* PITA_NLDRIVERIMPL_H */
