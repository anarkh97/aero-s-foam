#ifndef ROM_DISTRELEMENTSAMPLINGDRIVER_H
#define ROM_DISTRELEMENTSAMPLINGDRIVER_H

#include "DriverInterface.h"

class Communicator;

namespace Rom {

class DistrElementSamplingDriver : public DriverInterface {
public:
  virtual void solve();
  
  DistrElementSamplingDriver(Domain *, Communicator *);

private:
  Domain *domain_;
  Communicator *comm_;
};

} /* end namespace Rom */

Rom::DriverInterface *distrElementSamplingDriverNew(Domain *);

#endif /* ROM_DISTRELEMENTSAMPLINGDRIVER_H */
