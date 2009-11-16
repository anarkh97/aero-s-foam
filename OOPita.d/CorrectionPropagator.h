#ifndef PITA_CORRECTIONPROPAGATOR_H
#define PITA_CORRECTIONPROPAGATOR_H

#include "Fwk.h"
#include "NamedTask.h"

namespace Pita {

template <typename S> class SharedState;

template <typename S>
class CorrectionPropagator : public NamedTask {
public:
  EXPORT_PTRINTERFACE_TYPES(CorrectionPropagator);
  typedef Fwk::GenManagerInterface<CorrectionPropagator<S> *, String> Manager;
  
  // Inputs  
  const SharedState<S> * jump() const { return jump_.ptr(); }
  const SharedState<S> * correction() const { return correction_.ptr(); }
 
  virtual void jumpIs(const SharedState<S> * as) { setJump(as); }
  virtual void correctionIs(const SharedState<S> * sc) { setCorrection(sc); }

  // Outputs
  SharedState<S> * nextCorrection() const { return nextCorrection_.ptr(); }
 
  virtual void nextCorrectionIs(SharedState<S> * nc) { setNextCorrection(nc); }
  
protected:
  explicit CorrectionPropagator(const String & name) :
    NamedTask(name)
  {}
  
  void setCorrection(const SharedState<S> * c) { correction_ = c; }
  void setJump(const SharedState<S> * j) { jump_ = j; }
  void setNextCorrection(SharedState<S> * nc) { nextCorrection_ = nc; }
  
private:
  Fwk::Ptr<const SharedState<S> > jump_;
  Fwk::Ptr<const SharedState<S> > correction_;
  Fwk::Ptr<SharedState<S> > nextCorrection_;

  DISALLOW_COPY_AND_ASSIGN(CorrectionPropagator);
};
  
} /* end namespace Pita */

#endif /* PITA_CORRECTIONPROPAGATOR_H */
