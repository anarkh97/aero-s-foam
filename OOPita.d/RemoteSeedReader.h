#ifndef PITA_REMOTESEEDREADER_H
#define PITA_REMOTESEEDREADER_H

#include "Fwk.h"
#include "Types.h"
#include "HalfSliceTypes.h"

#include "Seed.h"

#include "SimpleBuffer.h"
class Communicator;

namespace Pita {

class RemoteSeedReader : public Seed::NotifieeConst {
public:
  EXPORT_PTRINTERFACE_TYPES(RemoteSeedReader);

  enum Status {
    READY = 0,
    SCHEDULED,
    BUSY
  };
  
  template <typename Key> class Manager;

  virtual void onState(); // Overriden

  CpuRank targetCpu() const { return targetCpu_; }
  Status status() const { return status_; }

  void targetCpuIs(CpuRank c) { targetCpu_ = c; }
  void statusIs(Status s);

protected:
  template <typename Key> friend class Manager;

  explicit RemoteSeedReader(Communicator * timeComm);

  void setStatus(Status s) { status_ = s; }

private:
  CpuRank targetCpu_;
  Status status_;

  Communicator * timeComm_;
  SimpleBuffer<double> stateBuffer_;
};


template <typename Key>
class RemoteSeedReader::Manager : public Fwk::PtrInterface<Manager<Key> >, private Fwk::GenManagerImpl<RemoteSeedReader, Key> {
public:
  EXPORT_PTRINTERFACE_TYPES(Manager<Key>);

  RemoteSeedReader * instance(const Key & key) const { return Fwk::GenManagerImpl<RemoteSeedReader, Key>::instance(key); } 
  size_t instanceCount() const { return Fwk::GenManagerImpl<RemoteSeedReader, Key>::instanceCount(); }

  RemoteSeedReader * instanceNew(const Key & key) { return Fwk::GenManagerImpl<RemoteSeedReader, Key>::instanceNew(key); }
  void instanceDel(const Key & key) { Fwk::GenManagerImpl<RemoteSeedReader, Key>::instanceDel(key); }

  static Ptr New(Communicator * comm) {
    return new Manager(comm);
  } 

protected:
  explicit Manager(Communicator * comm) :
    timeComm_(comm)
  {}

  virtual RemoteSeedReader * createNewInstance(const Key & key);

private:
  Communicator * timeComm_;
};

/* Template fix */
template <typename Key>
RemoteSeedReader *
RemoteSeedReader::Manager<Key>::createNewInstance(const Key & key) {
  return new RemoteSeedReader(timeComm_);;
}

} // end namespace Pita

#endif /* PITA_REMOTESEEDREADER_H */
