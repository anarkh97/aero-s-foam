#ifndef PITA_REMOTESTATEMPIIMPL_H
#define PITA_REMOTESTATEMPIIMPL_H

#include "RemoteState.h"
#include "Seed.h"

#include "SimpleBuffer.h"
#include <map>

class Communicator;

namespace Pita { namespace RemoteState {

class MpiManager;

/* Reader/Writer for Seed */

class MpiSeedReader : public Reader<DynamState> {
public:
  EXPORT_PTRINTERFACE_TYPES(MpiSeedReader);

  virtual void statusIs(Status s); // overriden

  Communicator * communicator() const { return communicator_; }

protected:
  MpiSeedReader(CpuRank targetCpu, Communicator * comm, SimpleBuffer<double> * buffer);

  const SimpleBuffer<double> & stateBuffer() const { return *stateBuffer_; }
  SimpleBuffer<double> & stateBuffer() { return *stateBuffer_; }

  friend class MpiManager;

private:
  Communicator * communicator_;
  SimpleBuffer<double> * stateBuffer_;
};

class MpiSeedWriter : public Writer<DynamState> {
public:
  EXPORT_PTRINTERFACE_TYPES(MpiSeedWriter);

  virtual void statusIs(Status s); // overriden

  Communicator * communicator() const { return communicator_; }

  size_t vectorSize() const { return vectorSize_; }
  void vectorSizeIs(size_t vSize) { vectorSize_ = vSize; }

protected:
  MpiSeedWriter(CpuRank originCpu, Communicator * comm, size_t vSize, SimpleBuffer<double> * buffer);

  const SimpleBuffer<double> & stateBuffer() const { return *stateBuffer_; }
  SimpleBuffer<double> & stateBuffer() { return *stateBuffer_; }

  friend class MpiManager;

private:
  Communicator * communicator_;
  size_t vectorSize_;
  SimpleBuffer<double> * stateBuffer_;
};

/* Reader/Writer for ReducedSeed */

class MpiReducedSeedReader : public Reader<Vector> {
public:
  EXPORT_PTRINTERFACE_TYPES(MpiReducedSeedReader);

  virtual void statusIs(Status s); // overriden

  Communicator * communicator() const { return communicator_; }

protected:
  MpiReducedSeedReader(CpuRank targetCpu, Communicator * comm, SimpleBuffer<double> * buffer);
  
  const SimpleBuffer<double> & stateBuffer() const { return *stateBuffer_; }
  SimpleBuffer<double> & stateBuffer() { return *stateBuffer_; }

  friend class MpiManager;

private:
  Communicator * communicator_;
  SimpleBuffer<double> * stateBuffer_;
};

class MpiReducedSeedWriter : public Writer<Vector> {
public:
  EXPORT_PTRINTERFACE_TYPES(MpiReducedSeedWriter);

  virtual void statusIs(Status s); // overriden

  Communicator * communicator() const { return communicator_; }

  size_t reducedStateSize() const { return reducedStateSize_; }
  void reducedStateSizeIs(size_t rSize) { reducedStateSize_ = rSize; }

protected:
  MpiReducedSeedWriter(CpuRank originCpu, Communicator * comm, size_t rSize, SimpleBuffer<double> * buffer);

  const SimpleBuffer<double> & stateBuffer() const { return *stateBuffer_; }
  SimpleBuffer<double> & stateBuffer() { return *stateBuffer_; }

  friend class MpiManager;

private:
  Communicator * communicator_;
  size_t reducedStateSize_;
  SimpleBuffer<double> * stateBuffer_;
};

/* MpiManager */

class MpiManager : public Manager {
public:
  /* Overriden */
  virtual MpiSeedReader * readerNew(const Seed * origin, CpuRank targetCpu);
  virtual MpiReducedSeedReader * readerNew(const ReducedSeed * origin, CpuRank targetCpu);

  virtual MpiSeedWriter * writerNew(Seed * target, CpuRank originCpu);
  virtual MpiReducedSeedWriter * writerNew(ReducedSeed * target, CpuRank originCpu);

  /* Added */
  Communicator * communicator() const { return communicator_; }

  size_t vectorSize() const { return vectorSize_; }
  size_t reducedStateSize() const { return reducedStateSize_; } // Assumed to be smaller than 2 * vectorSize

  void vectorSizeIs(size_t vSize);
  void reducedStateSizeIs(size_t rSize);

  static Ptr New(Communicator * comm) {
    return new MpiManager(comm);
  }

protected:
  explicit MpiManager(Communicator * comm);

private:
  Communicator * communicator_;
  
  size_t vectorSize_, reducedStateSize_;
  SimpleBuffer<double> sharedBuffer_;

  typedef std::multimap<const Seed *, MpiSeedReader::Ptr>               SeedReaderSet;
  typedef std::multimap<const Seed *, MpiSeedWriter::Ptr>               SeedWriterSet;
  typedef std::multimap<const ReducedSeed *, MpiReducedSeedReader::Ptr> ReducedSeedReaderSet;
  typedef std::multimap<const ReducedSeed *, MpiReducedSeedWriter::Ptr> ReducedSeedWriterSet;
  
  SeedReaderSet        seedReader_;
  SeedWriterSet        seedWriter_;
  ReducedSeedReaderSet reducedSeedReader_;
  ReducedSeedWriterSet reducedSeedWriter_;
};

} /* end namespace RemoteState */ } /* end namespace Pita */

#endif /* PITA_REMOTESTATEMPIIMPL_H */
