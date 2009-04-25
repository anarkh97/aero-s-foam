#ifndef PITA_NLDRIVER_H
#define PITA_NLDRIVER_H

#include "Fwk.h"

namespace Pita {

class PitaNonLinDynamic;
  
class NlDriver : public Fwk::PtrInterface<NlDriver> {
public:
  typedef Fwk::Ptr<NlDriver> Ptr;
  typedef Fwk::Ptr<const NlDriver> PtrConst;

  PitaNonLinDynamic * probDesc() const { return probDesc_; }

  // Main routine
  virtual void solve() = 0;

protected:
  explicit NlDriver(PitaNonLinDynamic * problemDescriptor) :
    probDesc_(problemDescriptor) {}

private:
  PitaNonLinDynamic * probDesc_;
};
  
} // end namespace Pita

extern Pita::NlDriver::Ptr nlPitaDriverNew(Pita::PitaNonLinDynamic * problemDescriptor);


#endif /* PITA_NLDRIVER_H */
