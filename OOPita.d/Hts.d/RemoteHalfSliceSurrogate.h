#ifndef PITA_REMOTEHALFSLICESURROGATE_H
#define PITA_REMOTEHALFSLICESURROGATE_H

#include "Fwk.h"
#include "Types.h"

#include "ScheduledRemoteSeedWriter.h"

#include "../Activity.h"
#include "SliceMapping.h"

namespace Pita {

// TODO HACK
using namespace Hts;

class RemoteHalfSliceSurrogate : public ScheduledRemoteSeedWriter::Scheduler {
public:
  EXPORT_PTRINTERFACE_TYPES(RemoteHalfSliceSurrogate);

  IterationRank iteration() const { return iteration_; }
  void iterationIs(IterationRank i);

  static RemoteHalfSliceSurrogate::Ptr New(SliceMapping * mapping) {
    return new RemoteHalfSliceSurrogate(mapping);
  }

protected:
  virtual void instanceNew(const Hts::SeedId & key, ScheduledRemoteSeedWriter *); // overriden
  virtual void instanceDel(const Hts::SeedId & key); // overriden

  explicit RemoteHalfSliceSurrogate(SliceMapping * mapping);

private:
  SliceMapping::Ptr mapping_;
  IterationRank iteration_;

  typedef ScheduledRemoteSeedWriter::ReceiveReactor ReceiveReactor;
  typedef std::map<Hts::SeedId, ReceiveReactor::Ptr> RecvRectorContainer;
  
  RecvRectorContainer evenReactor_;
  RecvRectorContainer oddReactor_;

  RecvRectorContainer & getReactorContainer(const Hts::SeedId & key);
  const RecvRectorContainer & getReactorContainer(const Hts::SeedId & key) const;
};


} // end namespace Pita

#endif /* PITA_REMOTEHALFSLICESURROGATE_H */
