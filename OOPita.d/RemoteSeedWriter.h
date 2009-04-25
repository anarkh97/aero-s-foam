#ifndef PITA_REMOTESEEDWRITER_H
#define PITA_REMOTESEEDWRITER_H

#include "Fwk.h"
#include "Types.h"
#include "HalfSliceTypes.h"

#include "Seed.h"

#include "SimpleBuffer.h"
class Communicator;

namespace Pita {

class RemoteSeedWriter : public Fwk::PtrInterface<RemoteSeedWriter> {
public:
  EXPORT_PTRINTERFACE_TYPES(RemoteSeedWriter);

  enum Status {
    READY = 0,
    BUSY
  };

  template <typename Key> class Manager;

  CpuRank originCpu() const { return originCpu_; } 
  Seed * target() const { return target_.ptr(); }
  size_t vectorSize() const { return vectorSize_; }
  Status status() const { return status_; }

  void originCpuIs(CpuRank oc) { setOriginCpu(oc); }
  void targetIs(Seed * t) { setTarget(t); }
  void vectorSizeIs(size_t vs) { setVectorSize(vs); }
  void statusIs(Status s);

protected:
  template <typename Key> friend class Manager;
  
  RemoteSeedWriter(Communicator * timeComm);

  void setOriginCpu(CpuRank oc) { originCpu_ = oc; }
  void setVectorSize(size_t vs) { vectorSize_ = vs; }
  void setStatus(Status s) { status_ = s; }
  void setTarget(Seed * t) { target_ = t; }

private:
  CpuRank originCpu_;
  Seed::Ptr target_;
  size_t vectorSize_;
  Status status_;

  Communicator * timeComm_;
  SimpleBuffer<double> stateBuffer_;
};

template <typename Key>
class RemoteSeedWriter::Manager : public Fwk::PtrInterface<RemoteSeedWriter::Manager<Key> >, private Fwk::GenManagerImpl<RemoteSeedWriter, Key> {
public:
  EXPORT_PTRINTERFACE_TYPES(Manager<Key>);
  typedef Fwk::GenManagerImpl<RemoteSeedWriter, Key> Impl;

  RemoteSeedWriter * instance(const Key & key) const { return Impl::instance(key); }
  size_t instanceCount() const { return Impl::instanceCount(); }
  
  RemoteSeedWriter * instanceNew(const Key & key) { return Impl::instanceNew(key); }
  void instanceDel(const Key & key) { Impl::instanceDel(key); } 

  static Ptr New(Communicator * timeComm) { 
    return new Manager(timeComm);
  }

protected:
  explicit Manager(Communicator * timeComm) :
    timeComm_(timeComm)
  {}

  virtual RemoteSeedWriter * createNewInstance(const Key & key);

private:
  Communicator * timeComm_;
};
 
template <typename Key>
RemoteSeedWriter *
RemoteSeedWriter::Manager<Key>::createNewInstance(const Key & key) {
  return new RemoteSeedWriter(timeComm_);
}

} // end namespace Pita

#endif
