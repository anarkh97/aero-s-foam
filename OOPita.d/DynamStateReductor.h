#ifndef PITA_DYNAMSTATEREDUCTOR_H
#define PITA_DYNAMSTATEREDUCTOR_H

#include "Fwk.h"

#include "DynamState.h"

#include <Math.d/Vector.h>

namespace Pita {

class DynamStateReductor : public Fwk::PtrInterface<DynamStateReductor> {
public:
  EXPORT_PTRINTERFACE_TYPES(DynamStateReductor);
  typedef Fwk::GenManagerInterface<DynamStateReductor*, String> Manager;

  size_t vectorSize() const { return vectorSize_; }

  const DynamState & initialState() const { return initialState_; } 
  virtual void initialStateIs(const DynamState & is) = 0;
  
  size_t reducedBasisSize() const { return reducedBasisComponents_.size(); }
  const Vector & reducedBasisComponents() const { return reducedBasisComponents_; }

protected:
  DynamStateReductor() :
    vectorSize_(0),
    initialState_(),
    reducedBasisComponents_()
  {}

  void setVectorSize(size_t vectorSize) { vectorSize_ = vectorSize; }
  void setInitialState(const DynamState & initialState) { initialState_ = initialState; }
  void setReducedBasisSize(size_t bSize) { reducedBasisComponents_.initialize(bSize); }

  Vector & getReducedBasisComponents() { return reducedBasisComponents_; } // Direct access to member for efficiency

private:
  size_t vectorSize_;
  DynamState initialState_;
  Vector reducedBasisComponents_;

  DISALLOW_COPY_AND_ASSIGN(DynamStateReductor);  
};

} // end namespace Pita

#endif /* PITA_DYNAMSTATEREDUCTOR_H */
