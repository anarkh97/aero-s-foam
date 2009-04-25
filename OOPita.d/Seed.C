#include "Seed.h"

namespace Pita {

Seed::Seed(const Fwk::String & name) :
  Fwk::NamedInterface(name),
  status_(Seed::INACTIVE)
{}

void
Seed::stateIs(const DynamState & s) {
  setState(s);
  notifierDelegate_.lastNotificationIs(&NotifieeConst::onState);
}

void
Seed::statusIs(Seed::Status s) {
  setStatus(s);
  notifierDelegate_.lastNotificationIs(&NotifieeConst::onStatus);
}

Seed::Manager::Manager() {
  // Nothing to do
}

Seed *
Seed::Manager::createNewInstance(const String & name) {
  return new Seed(name); 
}

} // end namespace Pita
