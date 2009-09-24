#ifndef PITA_HALFSLICENETWORK_H
#define PITA_HALFSLICENETWORK_H

#include <map>

#include "Fwk.h"

#include "../Activity.h"

#include "SliceMapping.h"
#include "HalfSliceNetworkTopology.h"

#include "../Seed.h"
#include "ScheduledRemoteSeedReader.h"
#include "ScheduledRemoteSeedWriter.h"

#include "HalfTimeSlice.h"
#include "FullTimeSlice.h"

#include "BasisCollector.h"

#include "HalfSliceSchedule.h"

#include "../SeedInitializer.h"

#include "RemoteHalfSliceSurrogate.h"

namespace Pita {

using namespace Hts;

class HalfSliceNetwork : public Fwk::PtrInterface<HalfSliceNetwork> {
public:
  EXPORT_PTRINTERFACE_TYPES(HalfSliceNetwork);

  HalfSliceNetwork(size_t vectorSize,
                   SliceMapping * mapping,
                   CpuRank localCpu,
                   Seed::Manager * seedMgr,
                   HalfTimeSlice::Manager * hsMgr,
                   FullTimeSliceHead::Manager * fshMgr,
                   FullTimeSliceTail::Manager * fstMgr,
                   ScheduledRemoteSeedReader::Manager<Hts::CommId> * sendMgr,
                   ScheduledRemoteSeedWriter::Manager * recvMgr,
                   RemoteHalfSliceSurrogate * surrogate, // HACK
                   HalfSliceSchedule * schedule,
                   SeedInitializer * primalSeedInit, 
                   SeedInitializer * dualSeedInit); 

protected:
  class SeedInitializationReactor : public Activity::Notifiee {
  public:
    typedef Fwk::Ptr<SeedInitializationReactor> Ptr;
    typedef Fwk::Ptr<const SeedInitializationReactor> PtrConst;

    virtual void onStatus(); // overriden

    HalfSliceNetwork * parent() const { return parent_; }

    SeedInitializationReactor(Activity * notifier, HalfSliceNetwork * parent);

  private:
    HalfSliceNetwork * parent_;
  };
  friend class SeedInitializationReactor;

  class NewIterationReactor : public Activity::Notifiee {
  public:
    typedef Fwk::Ptr<NewIterationReactor> Ptr;
    typedef Fwk::Ptr<const NewIterationReactor> PtrConst;

    virtual void onStatus(); // overriden

    HalfSliceNetwork * parent() const { return parent_; }

    NewIterationReactor(Activity * notifier, HalfSliceNetwork * parent);
  private:
    HalfSliceNetwork * parent_;
  };
  friend class NewIterationReactor;

  //friend class SeedIterator;

private:
  size_t vectorSize_;
  SliceMapping::Ptr mapping_;
 
  CpuRank localCpu_;

  Seed::Manager::Ptr seedMgr_;
  HalfTimeSlice::Manager::Ptr hsMgr_;
  FullTimeSliceHead::Manager::Ptr fshMgr_;
  FullTimeSliceTail::Manager::Ptr fstMgr_;
  
  //HalfSliceBasisCollector::Ptr collector_;

  ScheduledRemoteSeedReader::Manager<Hts::CommId>::Ptr sendMgr_;
  ScheduledRemoteSeedWriter::Manager::Ptr recvMgr_;
 
  RemoteHalfSliceSurrogate::Ptr surrogate_; // HACK

  HalfSliceSchedule::Ptr schedule_;
  HalfSliceNetworkTopology::PtrConst topology_;

  SeedInitializer::Ptr primalSeedInitializer_;
  SeedInitializer::Ptr dualSeedInitializer_;

  std::map<SliceRank, Seed *> primalInitialSeed_;
  std::map<SliceRank, Seed *> dualInitialSeed_;
  SeedInitializationReactor::Ptr seedInitializationReactor_;

  NewIterationReactor::Ptr newIterationReactor_;

  /*typedef std::map<HalfSliceRank, Seed::Ptr> SeedContainer;
  SeedContainer mainSeed_;
  SeedContainer leftPropagatedSeed_;
  SeedContainer rightPropagatedSeed_;

  SeedContainer & getSeedContainer(Hts::SeedType);*/
  void init();

  Seed * getSeed(const Hts::SeedId & id);
  PhaseRank getSeedSyncPhase(Hts::SliceType originSliceType, Hts::SliceType targetSliceType, HalfSliceRank remoteSliceRank) const;

  DISALLOW_COPY_AND_ASSIGN(HalfSliceNetwork);
};

/*class HalfSliceNetwork::Factory : public Fwk::PtrInterface<HalfSliceNetwork::Factory> {
public:
  typedef Fwk::Ptr<HalfSliceNetwork::Factory> Ptr;
  typedef Fwk::Ptr<const HalfSliceNetwork::Factory> PtrConst;
};*/

} // end namespace Pita

#endif /* PITA_HALFSLICENETWORK_H */
