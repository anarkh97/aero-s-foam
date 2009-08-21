#ifndef PITA_REDUCEDLINEARDRIVERIMPL_H
#define PITA_REDUCEDLINEARDRIVERIMPL_H

#include "LinearDriver.h"

#include "DynamState.h"

class Communicator;

namespace Pita {

class ReducedLinearDriverImpl : public LinearDriver {
public:
  EXPORT_PTRINTERFACE_TYPES(ReducedLinearDriverImpl);

  virtual void solve();

  SingleDomainDynamic<double> * probDesc() const { return probDesc_; }

  static ReducedLinearDriverImpl::Ptr New(SingleDomainDynamic<double> * pbDesc) {
    return new ReducedLinearDriverImpl(pbDesc);
  }

protected:
  explicit ReducedLinearDriverImpl(SingleDomainDynamic<double> * pbDesc);

private:
  SingleDomainDynamic<double> * probDesc_;

  /* Initial condition */
  size_t vectorSize_;
  DynamState initialSeed_;
 
  /* Time parameters */ 
  Seconds fineTimeStep_;
  TimeStepCount halfSliceRatio_;
  TimeStepCount sliceRatio_;
  Seconds coarseTimeStep_;
  Seconds initialTime_;
  Seconds finalTime_;
  
  /* Time domain decomposition */
  HalfSliceCount numSlices_;
  FullSliceCount fullTimeSlices_;

  /* CPUs & Network */ 
  Communicator * timeCom_;
  CpuRank myCpu_;
  CpuCount numCpus_;
  bool remoteCoarse_;
  HalfSliceCount maxActive_;
 
  /* Other paramters */ 
  IterationRank lastIteration_;
  double projectorTolerance_;
  bool noForce_;
};

Pita::LinearDriver::Ptr linearPitaDriverNew(SingleDomainDynamic<double> * pbDesc);

} // end namespace Pita

#endif /* PITA_REDUCEDLINEARDRIVERIMPL_H */
