#ifndef PITA_REMOTESHAREDSTATE_H
#define PITA_REMOTESHAREDSTATE_H

#include "Fwk.h"
#include "Types.h"

#include "SharedState.h"

#include "SimpleBuffer.h"
class Communicator;

namespace Pita {

/* Reader */

template <typename S>
class RemoteSharedStateReader : public SharedState<S>::NotifieeConst {
public:
  EXPORT_PTRINTERFACE_TYPES(RemoteSharedStateReader);

  enum Status {
    READY = 0,
    SCHEDULED,
    BUSY
  };
  
  template <typename Key> class Manager;

  virtual void onState(); // Overriden

  CpuRank targetCpu() const { return targetCpu_; }
  Status status() const { return status_; }

  virtual void targetCpuIs(CpuRank c) { targetCpu_ = c; }
  virtual void statusIs(Status s); // No default implementation, create specialization for each type S

protected:
  template <typename Key> friend class Manager;

  explicit RemoteSharedStateReader(Communicator * timeComm);

  void setStatus(Status s) { status_ = s; }

private:
  CpuRank targetCpu_;
  Status status_;

  Communicator * timeComm_;
  SimpleBuffer<double> stateBuffer_;
};


template <typename S>
template <typename Key>
class RemoteSharedStateReader<S>::Manager : public Fwk::PtrInterface<Manager<Key> >, private Fwk::GenManagerImpl<RemoteSharedStateReader<S>, Key> {
public:
  EXPORT_PTRINTERFACE_TYPES(Manager);

  RemoteSharedStateReader<S> * instance(const Key & key) const { return Fwk::GenManagerImpl<RemoteSharedStateReader<S>, Key>::instance(key); } 
  size_t instanceCount() const { return Fwk::GenManagerImpl<RemoteSharedStateReader<S>, Key>::instanceCount(); }

  RemoteSharedStateReader<S> * instanceNew(const Key & key) { return Fwk::GenManagerImpl<RemoteSharedStateReader<S>, Key>::instanceNew(key); }
  void instanceDel(const Key & key) { Fwk::GenManagerImpl<RemoteSharedStateReader<S>, Key>::instanceDel(key); }

  static Ptr New(Communicator * comm) {
    return new Manager(comm);
  } 

protected:
  explicit Manager(Communicator * comm) :
    timeComm_(comm)
  {}

  virtual RemoteSharedStateReader<S> * createNewInstance(const Key & key);

private:
  Communicator * timeComm_;
};

template <typename S>
template <typename Key>
RemoteSharedStateReader<S> *
RemoteSharedStateReader<S>::Manager<Key>::createNewInstance(const Key & key) {
  return new RemoteSharedStateReader<S>(timeComm_);;
}

template <typename S>
RemoteSharedStateReader<S>::RemoteSharedStateReader(Communicator * comm) :
  SharedState<S>::NotifieeConst(NULL),
  targetCpu_(),
  status_(READY),
  timeComm_(comm),
  stateBuffer_()
{}

template <typename S>
void
RemoteSharedStateReader<S>::onState() {
  this->statusIs(BUSY);
}


/* Writer */

template <typename S>
class RemoteSharedStateWriter : public Fwk::PtrInterface<RemoteSharedStateWriter<S> > {
public:
  EXPORT_PTRINTERFACE_TYPES(RemoteSharedStateWriter);

  enum Status {
    READY = 0,
    BUSY
  };

  template <typename Key> class Manager;

  CpuRank originCpu() const { return originCpu_; } 
  SharedState<S> * target() const { return target_.ptr(); }
  size_t vectorSize() const { return vectorSize_; }
  Status status() const { return status_; }

  virtual void originCpuIs(CpuRank oc) { setOriginCpu(oc); }
  virtual void targetIs(SharedState<S> * t) { setTarget(t); }
  virtual void vectorSizeIs(size_t vs) { setVectorSize(vs); }
  virtual void statusIs(Status s); // No generic implementation, provide specialization for each type S

protected:
  template <typename Key> friend class Manager;
  
  RemoteSharedStateWriter(Communicator * timeComm);

  void setOriginCpu(CpuRank oc) { originCpu_ = oc; }
  void setVectorSize(size_t vs) { vectorSize_ = vs; }
  void setStatus(Status s) { status_ = s; }
  void setTarget(SharedState<S> * t) { target_ = t; }

private:
  CpuRank originCpu_;
  typename SharedState<S>::Ptr target_;
  size_t vectorSize_;
  Status status_;

  Communicator * timeComm_;
  SimpleBuffer<double> stateBuffer_;
};

template <typename S>
RemoteSharedStateWriter<S>::RemoteSharedStateWriter(Communicator * timeComm) :
    originCpu_(),
    target_(NULL),
    vectorSize_(0),
    status_(READY),
    timeComm_(timeComm),
    stateBuffer_()
{}

/*
void
RemoteSharedState::statusIs(Status s) {
  if (s != status()) {
    setStatus(s);
    if (s == BUSY) {
      if (originCpu().value() != timeComm_->myID()) {
        log() << "RemoteWriter for " << (target() ? target()->name() : "None") << " from Cpu " << originCpu()  << " is busy\n";
        size_t dim = 2 * vectorSize_ + 1;
        stateBuffer_.sizeIs(dim);
        int messageTag = originCpu().value(); // TODO messageTag
        timeComm_->recFrom(messageTag, stateBuffer_.array(), dim); // Blocking MPI_RECV
        log() << "RemoteWriter for " << (target() ? target()->name() : "None") << " from Cpu " << originCpu()  << " is done\n";
        if (target()) {
          target()->statusIs(Seed::Status(static_cast<int>(stateBuffer_[dim - 1])));
          target()->stateIs(DynamState(vectorSize_, stateBuffer_.array()));
        }
      }
      setStatus(READY);
    }
  }
}
*/

template <typename S>
template <typename Key>
class RemoteSharedStateWriter<S>::Manager : public Fwk::PtrInterface<Manager<Key> >, private Fwk::GenManagerImpl<RemoteSharedStateWriter<S>, Key> {
public:
  EXPORT_PTRINTERFACE_TYPES(Manager<Key>);
  typedef Fwk::GenManagerImpl<RemoteSharedStateWriter<S>, Key> Impl;

  RemoteSharedStateWriter<S> * instance(const Key & key) const { return Impl::instance(key); }
  size_t instanceCount() const { return Impl::instanceCount(); }
  
  RemoteSharedStateWriter<S> * instanceNew(const Key & key) { return Impl::instanceNew(key); }
  void instanceDel(const Key & key) { Impl::instanceDel(key); } 

  static Ptr New(Communicator * timeComm) { 
    return new Manager(timeComm);
  }

protected:
  explicit Manager(Communicator * timeComm) :
    timeComm_(timeComm)
  {}

  virtual RemoteSharedStateWriter<S> * createNewInstance(const Key & key);

private:
  Communicator * timeComm_;
};
 
template <typename S>
template <typename Key>
RemoteSharedStateWriter<S> *
RemoteSharedStateWriter<S>::Manager<Key>::createNewInstance(const Key & key) {
  return new RemoteSharedStateWriter<S>(timeComm_);
}

} // end namespace Pita

#endif /* PITA_REMOTESHAREDSTATE_H */
