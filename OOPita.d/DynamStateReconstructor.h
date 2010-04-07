#ifndef PITA_DYNAMSTATERECONSTRUCTOR_H
#define PITA_DYNAMSTATERECONSTRUCTOR_H

#include "Fwk.h"

#include "DynamState.h"
#include "DynamPropagator.h"

#include <Math.d/Vector.h>

namespace Pita {

class DynamStateReconstructor : public DynamPropagatorTail {
public:
  EXPORT_PTRINTERFACE_TYPES(DynamStateReconstructor);
  typedef Fwk::GenManagerInterface<DynamStateReconstructor*, String> Manager;

  size_t reducedBasisSize() const { return reducedBasisSize_; }
  
  // Removed as an optimization
  // const Vector & reducedBasisComponents() const { return reducedBasisComponents_; }
  virtual void reducedBasisComponentsIs(const Vector & c) = 0;

protected:
  DynamStateReconstructor() :
    DynamPropagatorTail(0),
    reducedBasisSize_(0)
  {}

  void setReducedBasisSize(size_t size) { reducedBasisSize_ = size; }
 
private:
  size_t reducedBasisSize_;

  DISALLOW_COPY_AND_ASSIGN(DynamStateReconstructor);
};

} // end namespace Pita

#endif /* PITA_DYNAMSTATERECONSTRUCTOR_H */
