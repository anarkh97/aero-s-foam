#ifndef PITA_DYNAMSTATERECONSTRUCTOR_H
#define PITA_DYNAMSTATERECONSTRUCTOR_H

#include "Fwk.h"

#include "DynamState.h"
#include <Math.d/Vector.h>

#include "DynamStateBasis.h"

namespace Pita {

class DynamStateReconstructor : public Fwk::PtrInterface<DynamStateReconstructor> {
public:
  EXPORT_PTRINTERFACE_TYPES(DynamStateReconstructor);

  // Operator characteristics 
  size_t vectorSize() const { return reconstructionBasis_->vectorSize(); }
  size_t reducedBasisSize() const { return reconstructionBasis_->stateCount(); }
  
  const DynamStateBasis * reconstructionBasis() const { return reconstructionBasis_.ptr(); } 
  void reconstructionBasisIs(const DynamStateBasis * rb);

  // Input
  virtual void reducedBasisComponentsIs(const Vector & c);
  
  // Result
  const DynamState & finalState() const { return finalState_; }
  
  explicit DynamStateReconstructor(const DynamStateBasis * rb);

protected:
  void setFinalState(const DynamState & finalState) { finalState_ = finalState; }
 
private:
  DynamState finalState_;
  DynamStateBasis::PtrConst reconstructionBasis_;

  DISALLOW_COPY_AND_ASSIGN(DynamStateReconstructor);
};

} // end namespace Pita

#endif /* PITA_DYNAMSTATERECONSTRUCTOR_H */
