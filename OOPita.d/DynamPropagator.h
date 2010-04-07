#ifndef PITA_DYNAMPROPAGATOR_H
#define PITA_DYNAMPROPAGATOR_H

#include "Fwk.h"
#include "Types.h"
#include "DynamState.h"

namespace Pita {

class DynamPropagatorRoot : public Fwk::PtrInterface<DynamPropagatorRoot> {
public:
  EXPORT_PTRINTERFACE_TYPES(DynamPropagatorRoot);

  size_t vectorSize() const { return vectorSize_; }

protected:
  explicit DynamPropagatorRoot(size_t vectorSize = 0) :
    vectorSize_(vectorSize)
  {}

  void setVectorSize(size_t v) { vectorSize_ = v; }

private:
  size_t vectorSize_;

  DISALLOW_COPY_AND_ASSIGN(DynamPropagatorRoot);
};

class DynamPropagatorHead : virtual public DynamPropagatorRoot {
public:
  EXPORT_PTRINTERFACE_TYPES(DynamPropagatorHead);

  const DynamState & initialState() const { return initialState_; }
  virtual void initialStateIs(const DynamState & is) = 0;

protected:
  explicit DynamPropagatorHead(size_t vectorSize) :
    DynamPropagatorRoot(vectorSize),
    initialState_(vectorSize)
  {}

  void setInitialState(const DynamState & is) { initialState_ = is; }

private:
  DynamState initialState_;

  DISALLOW_COPY_AND_ASSIGN(DynamPropagatorHead);
};

class DynamPropagatorTail : virtual public DynamPropagatorRoot {
public:
  EXPORT_PTRINTERFACE_TYPES(DynamPropagatorTail);

  const DynamState & finalState() const { return finalState_; }

protected:
  explicit DynamPropagatorTail(size_t vectorSize) :
    DynamPropagatorRoot(vectorSize),
    finalState_(vectorSize)
  {}

  void setFinalState(const DynamState & fs) { finalState_ = fs; }

private:
  DynamState finalState_;

  DISALLOW_COPY_AND_ASSIGN(DynamPropagatorTail);
};

class DynamPropagator : public DynamPropagatorHead, public DynamPropagatorTail {
public:
  EXPORT_PTRINTERFACE_TYPES(DynamPropagator);

  // Overriden
  virtual void initialStateIs(const DynamState & is) { setInitialState(is); };

  class Notifiee : public Fwk::BaseMultiNotifiee<const DynamPropagator> {
  public:
    EXPORT_PTRINTERFACE_TYPES(Notifiee);

    virtual void onInitialState() {}
    virtual void onFinalState()   {}

  protected:
    Notifiee(const DynamPropagator * notifier) :
      Fwk::BaseMultiNotifiee<const DynamPropagator>(notifier)
    {}
  };

  // Notification interface
  virtual void lastNotifieeIs(Notifiee * notifiee) const { this->notifierDelegate().lastNotifieeIs(notifiee); }
  virtual void notifieeDel(Notifiee * notifiee) const { this->notifierDelegate().notifieeDel(notifiee); }

protected:
  explicit DynamPropagator(size_t vectorSize) :
    DynamPropagatorHead(vectorSize),
    DynamPropagatorTail(vectorSize)
  {
    setVectorSize(vectorSize);
  }

  void initialStateNotify() { this->notifierDelegate().lastNotificationIs(&Notifiee::onInitialState); }
  void finalStateNotify() { this->notifierDelegate().lastNotificationIs(&Notifiee::onFinalState); } 

  GenNotifierDelegate<Notifiee> & notifierDelegate() const { return const_cast<DynamPropagator *>(this)->notifierDelegate_; }

private:
  GenNotifierDelegate<Notifiee> notifierDelegate_;

  DISALLOW_COPY_AND_ASSIGN(DynamPropagator);
};

} // end namespace Pita

#endif /* PITA_DYNAMPROPAGATOR_H */
