#ifndef PITA_PROJECTORPROPAGATOR_H
#define PITA_PROJECTORPROPAGATOR_H

#include "Fwk.h"
#include "Types.h"
#include "DynamState.h"
#include "DynamStatePlainBasis.h"
#include "DynamPropagator.h"
#include "NlDynamOps.h"

namespace Pita {
  
class ProjectorPropagator : public DynamPropagator {
public:
  typedef Fwk::Ptr<ProjectorPropagator> Ptr;
  typedef Fwk::Ptr<const ProjectorPropagator> PtrConst;

  // Base attributes 
  virtual void initialStateIs(const DynamState & initialState);
  
  // Specialized attributes
  const GenVector<double> & referenceDisplacement() const { return refDisp_; }
  void referenceDisplacementIs(const GenVector<double> & disp, Seconds time);
  
  TimeStepCount timeStepCount() const { return timeStepCount_; }
  void timeStepCountInc();
 
  DynamStatePlainBasis::PtrConst projectionBasis() const { return projectionBasis_; }
  void projectionBasisIs(DynamStateBasis::PtrConst basis);
  void projectionBasisInc(DynamStateBasis::PtrConst lastBasis);
  
  DynamStatePlainBasis::PtrConst propagatedBasis() const { return propagatedBasis_; }
  
  // Instanciation
  static Ptr New(DynamPropagator * basePropagator, NlDynamOps * dynamOps) {
    return new ProjectorPropagator(basePropagator, dynamOps);
  }
  
  // Debug 
  void toString(OStream &) const; 
  
protected:
  ProjectorPropagator(DynamPropagator * basePropagator, NlDynamOps * dynamOps);
  
private:
  static const double tolerance_;
  
  DynamPropagator::Ptr basePropagator_;
  NlDynamOps::Ptr dynamOps_;
  GenVector<double> refDisp_;
  Seconds refTime_;
  TimeStepCount timeStepCount_;
  
  DynamStatePlainBasis::Ptr projectionBasis_;
  DynamStatePlainBasis::Ptr orthoBasis_;
  DynamStatePlainBasis::Ptr propagatedBasis_;
};

} // end namespace Pita

#endif /* PITA_PROJECTORPROPAGATOR_H */
