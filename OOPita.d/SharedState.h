#ifndef PITA_SHAREDSTATE_H
#define PITA_SHAREDSTATE_H

#include "Fwk.h"

namespace Pita {

template <typename S>
class SharedState : public Fwk::NamedInterface {
public:
  EXPORT_PTRINTERFACE_TYPES(SharedState);

  typedef S StateType;
  class Manager;
  
  enum Status {
    INACTIVE = 0,
    ACTIVE,
    CONVERGED,
    SPECIAL
  };

  const StateType & state() const { return state_; }
  Status status() const { return status_; }

  void stateIs(const StateType & s); 
  void statusIs(Status s);

  class NotifieeConst : public Fwk::BaseMultiNotifiee<const SharedState, NotifieeConst> {
  public:
    typedef Fwk::Ptr<NotifieeConst> Ptr;
    typedef Fwk::Ptr<const NotifieeConst> PtrConst;

    virtual void onState() {}
    virtual void onStatus() {}

  protected:
    explicit NotifieeConst(const SharedState * notifier = NULL) :
      Fwk::BaseMultiNotifiee<const SharedState, NotifieeConst>(notifier)
    {}
  };

  // Notifier implementation
  void lastNotifieeIs(NotifieeConst * notifiee) const { notifierDelegate().lastNotifieeIs(notifiee); }
  void notifieeDel(NotifieeConst * notifiee) const { notifierDelegate().notifieeDel(notifiee); }

protected:
  explicit SharedState(const Fwk::String & name);

  void setState(const StateType & s) { state_ = s; }
  void setStatus(Status s) { status_ = s; }
  
  GenNotifierDelegate<NotifieeConst> & notifierDelegate() const { return const_cast<SharedState *>(this)->notifierDelegate_; }
  
  friend class Manager;

private:
  StateType state_;
  Status status_;
  
  // Notifier implementation
  GenNotifierDelegate<NotifieeConst> notifierDelegate_;

  DISALLOW_COPY_AND_ASSIGN(SharedState);
};

template <typename S>
class SharedState<S>::Manager : public Fwk::GenNamedInterfaceManager<SharedState<S> > {
public:
  EXPORT_PTRINTERFACE_TYPES(Manager);

  static Ptr New() {
    return new Manager();
  }

protected:
  Manager();

  virtual SharedState<S> * createNewInstance(const String & name); // Overriden
};


template <typename S>
SharedState<S>::SharedState(const Fwk::String & name) :
  Fwk::NamedInterface(name),
  status_(SharedState<S>::INACTIVE)
{}

template <typename S>
void
SharedState<S>::stateIs(const SharedState<S>::StateType & s) {
  setState(s);
  notifierDelegate_.lastNotificationIs(&NotifieeConst::onState);
}

template <typename S>
void
SharedState<S>::statusIs(SharedState<S>::Status s) {
  setStatus(s);
  notifierDelegate_.lastNotificationIs(&NotifieeConst::onStatus);
}

template <typename S>
SharedState<S>::Manager::Manager() {
  // Nothing to do
}

template <typename S>
SharedState<S> *
SharedState<S>::Manager::createNewInstance(const String & name) {
  return new SharedState<S>(name); 
}

} // end namespace Pita

#endif /* PITA_SHAREDSTATE_H */
