#ifndef PITA_REMOTEPROPAGATOR_H
#define PITA_REMOTEPROPAGATOR_H

#include "DynamPropagator.h"
#include "SimpleBuffer.h"
#include <Comm.d/Communicator.h>

namespace Pita {

class RemoteDynamPropagator : public DynamPropagator {
public:
  EXPORT_PTRINTERFACE_TYPES(RemoteDynamPropagator);

  // Added
  Communicator * serverCommunicator() const { return serverCommunicator_; }
  CpuRank serverCpu() const { return serverCpu_; }
  
  // Overriden
  virtual void initialStateIs(const DynamState & is);

  static Ptr New(size_t vectorSize, Communicator * serverCommunicator, CpuRank serverCpu) {
    return new RemoteDynamPropagator(vectorSize, serverCommunicator, serverCpu);
  }

protected:
  RemoteDynamPropagator(size_t vectorSize, Communicator * serverCommunicator, CpuRank serverCpu);
  
private:
  Communicator * serverCommunicator_;
  CpuRank serverCpu_; 

  SimpleBuffer<double> sBuffer_;
};

} // end namespace Pita

#endif /* PITA_REMOTEPROPAGATOR_H */
