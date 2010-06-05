#ifndef PITA_DYNAMSTATERECONSTRUCTOR_H
#define PITA_DYNAMSTATERECONSTRUCTOR_H

#include "Fwk.h"

#include "DynamState.h"
#include <Math.d/Vector.h>

namespace Pita {

class DynamStateReconstructor : public Fwk::PtrInterface<DynamStateReconstructor> {
public:
  EXPORT_PTRINTERFACE_TYPES(DynamStateReconstructor);
  typedef Fwk::GenManagerInterface<DynamStateReconstructor*, String> Manager;

  size_t vectorSize() const { return vectorSize_; }
  const DynamState & finalState() const { return finalState_; }
  size_t reducedBasisSize() const { return reducedBasisSize_; }
  
  // Removed as an optimization
  // const Vector & reducedBasisComponents() const { return reducedBasisComponents_; }
  virtual void reducedBasisComponentsIs(const Vector & c) = 0;

protected:
  DynamStateReconstructor() :
    vectorSize_(0),
    finalState_(),
    reducedBasisSize_(0)
  {}

  void setVectorSize(size_t vectorSize) { vectorSize_ = vectorSize; }
  void setFinalState(const DynamState & finalState) { finalState_ = finalState; }
  void setReducedBasisSize(size_t size) { reducedBasisSize_ = size; }
 
private:
  size_t vectorSize_;
  DynamState finalState_;
  size_t reducedBasisSize_;

  DISALLOW_COPY_AND_ASSIGN(DynamStateReconstructor);
};

} // end namespace Pita

#endif /* PITA_DYNAMSTATERECONSTRUCTOR_H */
