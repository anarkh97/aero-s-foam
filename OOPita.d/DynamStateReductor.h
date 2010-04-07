#ifndef PITA_DYNAMSTATEREDUCTOR_H
#define PITA_DYNAMSTATEREDUCTOR_H

#include "Fwk.h"

#include "DynamState.h"
#include "DynamPropagator.h"

#include <Math.d/Vector.h>

namespace Pita {

class DynamStateReductor : public DynamPropagatorHead {
public:
  EXPORT_PTRINTERFACE_TYPES(DynamStateReductor);
  typedef Fwk::GenManagerInterface<DynamStateReductor*, String> Manager;
  
  // Overriden  
  virtual void initialStateIs(const DynamState & is) = 0;
  
  size_t reducedBasisSize() const { return reducedBasisComponents_.size(); }
  const Vector & reducedBasisComponents() const { return reducedBasisComponents_; }

protected:
  DynamStateReductor() :
    DynamPropagatorHead(0),
    reducedBasisComponents_()
  {}

  void setReducedBasisSize(size_t bSize) { reducedBasisComponents_.initialize(bSize); }

  Vector & getReducedBasisComponents() { return reducedBasisComponents_; } // Direct access to member for efficiency

private:
  Vector reducedBasisComponents_;

  DISALLOW_COPY_AND_ASSIGN(DynamStateReductor);  
};

} // end namespace Pita

#endif /* PITA_DYNAMSTATEREDUCTOR_H */
