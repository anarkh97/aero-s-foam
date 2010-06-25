#ifndef PITA_NLDRIVER_H
#define PITA_NLDRIVER_H

#include "Fwk.h"

namespace Pita {

class PitaNonLinDynamic;
  
class NlDriver : public Fwk::PtrInterface<NlDriver> {
public:
  EXPORT_PTRINTERFACE_TYPES(NlDriver);

  PitaNonLinDynamic * probDesc() const { return probDesc_; }

  // Main routine
  virtual void solve() = 0;

protected:
  explicit NlDriver(PitaNonLinDynamic * problemDescriptor) :
    probDesc_(problemDescriptor) {}

private:
  PitaNonLinDynamic * probDesc_;

  DISALLOW_COPY_AND_ASSIGN(NlDriver);
};
  
} // end namespace Pita

extern Pita::NlDriver::Ptr nlPitaDriverNew(Pita::PitaNonLinDynamic * problemDescriptor);


#endif /* PITA_NLDRIVER_H */
