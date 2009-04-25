#ifndef PITA_SEED_H
#define PITA_SEED_H

#include "Fwk.h"
#include "DynamState.h"
#include "Types.h"

namespace Pita {

class Seed : public Fwk::NamedInterface {
public:
  EXPORT_PTRINTERFACE_TYPES(Seed);

  class Manager;
  
  enum Status {
    INACTIVE = 0,
    ACTIVE,
    CONVERGED,
    SPECIAL
  };

  const DynamState & state() const { return state_; }
  Status status() const { return status_; }

  void stateIs(const DynamState & s); 
  void statusIs(Status s);

  class NotifieeConst : public Fwk::BaseMultiNotifiee<const Seed, NotifieeConst> {
  public:
    typedef Fwk::Ptr<NotifieeConst> Ptr;
    typedef Fwk::Ptr<const NotifieeConst> PtrConst;

    virtual void onState() {}
    virtual void onStatus() {}

  protected:
    explicit NotifieeConst(const Seed * notifier = NULL) :
      Fwk::BaseMultiNotifiee<const Seed, NotifieeConst>(notifier)
    {}
  };

  // Notifier implementation
  void lastNotifieeIs(NotifieeConst * notifiee) const { notifierDelegate().lastNotifieeIs(notifiee); }
  void notifieeDel(NotifieeConst * notifiee) const { notifierDelegate().notifieeDel(notifiee); }

protected:
  explicit Seed(const Fwk::String & name);

  void setState(const DynamState & s) { state_ = s; }
  void setStatus(Status s) { status_ = s; }
  
  GenNotifierDelegate<NotifieeConst> & notifierDelegate() const { return const_cast<Seed *>(this)->notifierDelegate_; }
  
  friend class Manager;

private:
  DynamState state_;
  Status status_;
  
  // Notifier implementation
  GenNotifierDelegate<NotifieeConst> notifierDelegate_;

  DISALLOW_COPY_AND_ASSIGN(Seed);
};


class Seed::Manager : public Fwk::GenNamedInterfaceManager<Seed> {
public:
  EXPORT_PTRINTERFACE_TYPES(Manager);

  static Ptr New() {
    return new Manager();
  }

protected:
  Manager();

  virtual Seed * createNewInstance(const String & name); // Overriden
};

} // end namespace Pita

#endif /* PITA_SEED_H */
