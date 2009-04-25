#ifndef PITA_HALFSLICEBASISCOLLECTORIMPL_H
#define PITA_HALFSLICEBASISCOLLECTORIMPL_H

#include "HalfSliceBasisCollector.h"

#include "stack"

namespace Pita {

class HalfSliceBasisCollectorImpl : public HalfSliceBasisCollector {
public: 
  EXPORT_PTRINTERFACE_TYPES(HalfSliceBasisCollectorImpl);
  
  // Overriden interface
  virtual IntegratorPropagator * source(const HalfSliceId & sliceId) const;
  virtual size_t sourceCount() const;
  virtual IntegratorPropagator * sourceNew(const HalfSliceId & sliceId, DynamTimeIntegrator * baseIntegrator);
  virtual void sourceDel(const HalfSliceId & sliceId);

  // Specialized interface 
  typedef std::pair<HalfSliceRank, DynamState> CollectedState;

  CollectedState firstForwardFinalState() const;
  void firstForwardFinalStateDel();
  CollectedState firstBackwardFinalState() const;
  void firstBackwardFinalStateDel();
  
  void finalStateIs(const HalfSliceId & sliceId, const DynamState & state); 

  Ptr static New() {
    return new HalfSliceBasisCollectorImpl();
  }

protected:
  HalfSliceBasisCollectorImpl();
 
  class PropagationReactor; 
  typedef std::map<HalfSliceId, Fwk::Ptr<PropagationReactor> > PropagationReactorContainer;
  
  const PropagationReactorContainer & propagationReactor() const {
    return propagationReactor_;
  }

  virtual PropagationReactor * propagationReactorNew(IntegratorPropagator * notifier,
                                                     const HalfSliceId & id);

private:
  PropagationReactorContainer propagationReactor_;
  
  typedef std::stack<CollectedState> StateContainer;
  StateContainer forwardFinalState_;
  StateContainer backwardFinalState_;
};

class HalfSliceBasisCollectorImpl::PropagationReactor : public DynamPropagator::Notifiee {
public:
  EXPORT_PTRINTERFACE_TYPES(PropagationReactor);
  
  PropagationReactor(DynamPropagator * notifier, const HalfSliceId & id, HalfSliceBasisCollectorImpl * parent);

  virtual void onFinalState(); // Overriden

  const HalfSliceId & sliceId() const { return sliceId_; }
  HalfSliceBasisCollectorImpl * parent() const { return parent_; }

private:
  HalfSliceId sliceId_;
  HalfSliceBasisCollectorImpl * parent_;
};

} // end namespace Pita

#endif /* PITA_HALFSLICEBASISCOLLECTORIMPL_H */
