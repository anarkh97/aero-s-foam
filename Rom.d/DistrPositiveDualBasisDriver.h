#ifndef ROM_DISTRPOSITIVEDUALBASISDRIVER_H
#define ROM_DISTRPOSITIVEDUALBASISDRIVER_H

#include "DriverInterface.h"

class Communicator;
class DistrInfo;

namespace Rom {

class DistrPositiveDualBasisDriver : public MultiDomainDynam, public DriverInterface {
public:
  virtual void solve();
  
  DistrPositiveDualBasisDriver(Domain *, Communicator *);

private:
  Domain *domain_;
  Communicator *comm_;
};

} /* end namespace Rom */

Rom::DriverInterface *distrPositiveDualBasisDriverNew(Domain *);

#endif /* ROM_DISTRPOSITIVEDUALBASISDRIVER_H */
