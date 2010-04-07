#ifndef PITA_NLDYNAMPOSTPROCESSOR_H
#define PITA_NLDYNAMPOSTPROCESSOR_H

#include "Fwk.h"
#include <OOPita.d/DynamPostProcessor.h>
#include <OOPita.d/NlDynamTimeIntegrator.h>

class GeomState;

namespace Pita {

class PitaNonLinDynamic;
  
class NlDynamPostProcessor : public DynamPostProcessor {
public:
  EXPORT_PTRINTERFACE_TYPES(NlDynamPostProcessor);

  virtual void statusIs(Status s);
  virtual void sliceRankIs(SliceRank r);
  
  void lastNlOutputIs(Seconds time, TimeStepCount step, const DynamState & state, const GeomState * geom, const GenVector<double> & force);

  static NlDynamPostProcessor::Ptr New(PitaNonLinDynamic * pbDesc) { 
    return new  NlDynamPostProcessor(pbDesc);
  }
  
  class IntegratorReactor : public DynamPostProcessor::IntegratorReactor {
  public:
    typedef Fwk::Ptr<IntegratorReactor> Ptr;
    typedef Fwk::Ptr<const IntegratorReactor> PtrConst;

    NlDynamPostProcessor::Ptr parent() const { return getParent(); }

    virtual void onInitialCondition();
    virtual void onCurrentCondition();
    
    NlDynamTimeIntegrator::Ptr nlNotifier() const { return nlNotifier_; }

    void nlNotifierIs(NlDynamTimeIntegrator::Ptr notifier) {
      notifierIs(notifier.ptr());
      nlNotifier_ = notifier;
    }
    
    friend class NlDynamPostProcessor;

  protected:
    IntegratorReactor(NlDynamTimeIntegrator * notifier, NlDynamPostProcessor * parent);

    virtual NlDynamPostProcessor * getParent() const { return parent_; }

  private:
    NlDynamTimeIntegrator::Ptr nlNotifier_;
    NlDynamPostProcessor * parent_;
  };

  IntegratorReactor::Ptr integratorReactorNew(DynamTimeIntegrator * notifier) {
    return getNewIntegratorReactor(notifier);
  }

  IntegratorReactor::Ptr integratorReactorNew(NlDynamTimeIntegrator * notifier);

protected:
  explicit NlDynamPostProcessor(PitaNonLinDynamic * pbDesc);

  ~NlDynamPostProcessor();
  
  virtual IntegratorReactor * getNewIntegratorReactor(const DynamTimeIntegrator * notifier);

  IntegratorReactor * getNewIntegratorReactor(NlDynamTimeIntegrator * notifier) {
    return new IntegratorReactor(notifier, this);
  }
  
private:
  PitaNonLinDynamic * probDesc_;
  GenVector<double> dummy_;
};
  
} // end namespace Pita

#endif /* PITA_NLDYNAMPOSTPROCESSOR_H */
