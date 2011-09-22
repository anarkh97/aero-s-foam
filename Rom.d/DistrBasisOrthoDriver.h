#ifndef ROM_DISTRBASISORTHODRIVER_H
#define ROM_DISTRBASISORTHODRIVER_H

#include "DriverInterface.h"

class Communicator;

namespace Rom {

class DistrBasisOrthoDriver : public DriverInterface {
public:
  virtual void solve();
  
  DistrBasisOrthoDriver(Domain *, Communicator *);

private:
  Domain *domain_;
  Communicator *comm_;
};

} /* end namespace Rom */

Rom::DriverInterface *distrBasisOrthoDriverNew(Domain *);

#endif /* ROM_DISTRBASISORTHODRIVER_H */
