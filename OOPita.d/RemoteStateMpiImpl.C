#include "RemoteStateMpiImpl.h"

#include "SimpleBuffer.h"
#include <Comm.d/Communicator.h>

namespace Pita { namespace RemoteState {

/* MpiSeedReader implementation */

MpiSeedReader::MpiSeedReader(CpuRank targetCpu, Communicator * comm, SimpleBuffer<double> * buffer) :
  Reader<DynamState>(targetCpu),
  communicator_(comm),
  stateBuffer_(buffer)
{}

void
MpiSeedReader::statusIs(Status s) {
  if (status() == s)
    return;

  if (s != INACTIVE && !isReady())
    throw Fwk::RangeException("in MpiSeedReader::statusIs -- Not ready");

  setStatus(s);

  if (status() != EXECUTING)
    return;

  if (communicator_->myID() != targetCpu().value() && targetCpu() != CpuRank(-1)) {
    DynamState state = origin()->state();

    size_t dim = 2 * state.vectorSize() + 2;
    bufferStateCopy(state, stateBuffer().array());

    stateBuffer()[dim - 2] = static_cast<double>(origin()->status()); 
    stateBuffer()[dim - 1] = static_cast<double>(origin()->iteration().value());

    int messageTag = communicator()->myID(); // TODO determine meaningful messageId
    //log() << "RemoteReader for " << origin()->name() << " to Cpu " << targetCpu()  << " is busy\n";
    //log() << "dim = " << dim << "\n";
    communicator()->sendTo(targetCpu().value(), messageTag, stateBuffer().array(), dim);
    communicator()->waitForAllReq(); // TODO delay waitForAllReq
    //log() << "RemoteReader for " << origin()->name() << " to Cpu " << targetCpu()  << " is done\n";
  }

  setStatus(READY);
}

/* MpiSeedWriter implementation */

MpiSeedWriter::MpiSeedWriter(CpuRank originCpu, Communicator * comm, size_t vSize, SimpleBuffer<double> * buffer) :
  Writer<DynamState>(originCpu),
  communicator_(comm),
  vectorSize_(vSize),
  stateBuffer_(buffer)
{}

void
MpiSeedWriter::statusIs(Status s) {
  if (status() == s)
    return;

  if (s != INACTIVE && !isReady())
    throw Fwk::RangeException("in MpiSeedWriter::statusIs -- Not ready");

  setStatus(s);
  
  if (status() != EXECUTING)
    return;

  if (originCpu().value() != communicator()->myID()) {
    //log() << "RemoteWriter for " << target()->name() << " from Cpu " << originCpu()  << " is busy\n";
    size_t dim = 2 * vectorSize() + 2;
    //log() << "dim = " << dim << "\n";
    int messageTag = originCpu().value(); // TODO messageTag
    communicator()->recFrom(messageTag, stateBuffer().array(), dim); // Blocking MPI_RECV
    //log() << "RemoteWriter for " << target()->name() << " from Cpu " << originCpu()  << " is done\n";
    target()->statusIs(Seed::Status(static_cast<int>(stateBuffer()[dim - 2])));
    target()->stateIs(DynamState(vectorSize(), stateBuffer().array()));
    target()->iterationIs(IterationRank(static_cast<int>(stateBuffer()[dim - 1])));
  }

  setStatus(READY);
}

/* MpiReducedSeedReader implementation */

MpiReducedSeedReader::MpiReducedSeedReader(CpuRank targetCpu, Communicator * comm, SimpleBuffer<double> * buffer) :
  Reader<Vector>(targetCpu),
  communicator_(comm),
  stateBuffer_(buffer)
{}

void
MpiReducedSeedReader::statusIs(Status s) {
  if (status() == s)
    return;

  if (s != INACTIVE && !isReady())
    throw Fwk::RangeException("in MpiReducedSeedReader::statusIs -- Not ready");

  setStatus(s);

  if (status() != EXECUTING)
    return;

  if (communicator_->myID() != targetCpu().value() && targetCpu() != CpuRank(-1)) {
    size_t dim = origin()->state().size() + 2;
    const double * stateBegin = origin()->state().data();
    std::copy(stateBegin, stateBegin + origin()->state().size(), stateBuffer().array());

    stateBuffer()[dim - 2] = static_cast<double>(origin()->status()); 
    stateBuffer()[dim - 1] = static_cast<double>(origin()->iteration().value());

    int messageTag = communicator()->myID(); // TODO determine meaningful messageId
    //log() << "RedRemoteReader for " << origin()->name() << " to Cpu " << targetCpu()  << " is busy\n";
    //log() << "dim = " << dim << "\n";
    communicator()->sendTo(targetCpu().value(), messageTag, stateBuffer().array(), dim);
    communicator()->waitForAllReq(); // TODO delay waitForAllReq
    //log() << "RedRemoteReader for " << origin()->name() << " to Cpu " << targetCpu()  << " is done\n";
  }

  setStatus(READY);
}

/* MpiReducedSeedWriter implementation */

MpiReducedSeedWriter::MpiReducedSeedWriter(CpuRank originCpu, Communicator * comm, size_t rSize, SimpleBuffer<double> * buffer) :
  Writer<Vector>(originCpu),
  communicator_(comm),
  reducedStateSize_(rSize),
  stateBuffer_(buffer)
{}

void
MpiReducedSeedWriter::statusIs(Status s) {
  if (status() == s)
    return;

  if (s != INACTIVE && !isReady())
    throw Fwk::RangeException("in MpiReducedSeedWriter::statusIs -- Not ready");

  setStatus(s);
  
  if (status() != EXECUTING)
    return;

  if (originCpu().value() != communicator()->myID()) {
    //log() << "RedRemoteWriter for " << target()->name() << " from Cpu " << originCpu()  << " is busy\n";
    size_t dim = reducedStateSize() + 2;
    //log() << "dim = " << dim << "\n";
    int messageTag = originCpu().value(); // TODO messageTag
    communicator()->recFrom(messageTag, stateBuffer().array(), dim); // Blocking MPI_RECV
    //log() << "RedRemoteWriter for " << target()->name() << " from Cpu " << originCpu()  << " is done\n";
    target()->statusIs(Seed::Status(static_cast<int>(stateBuffer()[dim - 2])));
    target()->stateIs(Vector(reducedStateSize(), stateBuffer().array()));
    target()->iterationIs(IterationRank(static_cast<int>(stateBuffer()[dim - 1])));
  }

  setStatus(READY);
}

/* MpiManager implementation */

MpiManager::MpiManager(Communicator * comm) :
  communicator_(comm),
  vectorSize_(0),
  reducedStateSize_(0),
  sharedBuffer_(0)
{}

MpiSeedReader *
MpiManager::readerNew(const Seed * origin, CpuRank targetCpu) {
  MpiSeedReader::Ptr reader = new MpiSeedReader(targetCpu, communicator(), &sharedBuffer_);
  reader->originIs(origin);
  seedReader_.insert(std::make_pair(origin, reader));

  return reader.ptr();
}

MpiSeedWriter *
MpiManager::writerNew(Seed * target, CpuRank originCpu) {
  MpiSeedWriter::Ptr writer = new MpiSeedWriter(originCpu, communicator(), vectorSize(), &sharedBuffer_);
  writer->targetIs(target);
  seedWriter_.insert(std::make_pair(const_cast<const Seed *>(target), writer));

  return writer.ptr();
}

MpiReducedSeedReader *
MpiManager::readerNew(const ReducedSeed * origin, CpuRank targetCpu) {
  MpiReducedSeedReader::Ptr reader = new MpiReducedSeedReader(targetCpu, communicator(), &sharedBuffer_);
  reader->originIs(origin);
  reducedSeedReader_.insert(std::make_pair(origin, reader));

  return reader.ptr();
}

MpiReducedSeedWriter *
MpiManager::writerNew(ReducedSeed * target, CpuRank originCpu) {
  MpiReducedSeedWriter::Ptr writer = new MpiReducedSeedWriter(originCpu, communicator(), reducedStateSize(), &sharedBuffer_);
  writer->targetIs(target);
  reducedSeedWriter_.insert(std::make_pair(const_cast<const ReducedSeed *>(target), writer));

  return writer.ptr();
}

void
MpiManager::vectorSizeIs(size_t vSize) {
  sharedBuffer_.sizeIs(2 * vSize + 2);

  SeedWriterSet::iterator it_end = seedWriter_.end();
  for (SeedWriterSet::iterator it = seedWriter_.begin(); it != it_end; ++it) {
    it->second->vectorSizeIs(vSize);
  }
  
  vectorSize_ = vSize;
}

void
MpiManager::reducedStateSizeIs(size_t rSize) {
  
  ReducedSeedWriterSet::iterator it_end = reducedSeedWriter_.end();
  for (ReducedSeedWriterSet::iterator it = reducedSeedWriter_.begin(); it != it_end; ++it) {
    it->second->reducedStateSizeIs(rSize);
  }
  
  reducedStateSize_ = rSize;
}

} /* end namespace RemoteState */ } /* end namespace Pita */
