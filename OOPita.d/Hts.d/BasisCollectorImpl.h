#ifndef PITA_HTS_BASISCOLLECTORIMPL_H
#define PITA_HTS_BASISCOLLECTORIMPL_H

#include "BasisCollector.h"

#include "stack"

namespace Pita { namespace Hts {

class BasisCollectorImpl : public BasisCollector {
public: 
  EXPORT_PTRINTERFACE_TYPES(BasisCollectorImpl);
  
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
    return new BasisCollectorImpl();
  }

protected:
  BasisCollectorImpl();
 
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

class BasisCollectorImpl::PropagationReactor : public DynamPropagator::Notifiee {
public:
  EXPORT_PTRINTERFACE_TYPES(PropagationReactor);
  
  PropagationReactor(DynamPropagator * notifier, const HalfSliceId & id, BasisCollectorImpl * parent);

  virtual void onFinalState(); // Overriden

  const HalfSliceId & sliceId() const { return sliceId_; }
  BasisCollectorImpl * parent() const { return parent_; }

private:
  HalfSliceId sliceId_;
  BasisCollectorImpl * parent_;
};

} /* end notifier Hts */ } /* end namespace Pita */

#endif /* PITA_HTS_BASISCOLLECTORIMPL_H */
