#ifndef ROM_DRIVERINTERFACE_H
#define ROM_DRIVERINTERFACE_H

#include "Paral.d/MDDynam.h"

class Domain;

namespace Rom {

class DriverInterface {
public:
  virtual void solve() = 0;
  virtual ~DriverInterface() {}

protected:
  DriverInterface() {}

private:
  // Disallow copy & assignment
  DriverInterface(const DriverInterface &);
  DriverInterface &operator=(const DriverInterface &);
};

} /* end namespace Rom */

// Concrete class instantiation
extern Rom::DriverInterface *basisOrthoDriverNew(Domain *);
extern Rom::DriverInterface *deimSamplingDriverNew(Domain *);
extern Rom::DriverInterface *udeimSamplingDriverNew(Domain *);
extern Rom::DriverInterface *distrBasisOrthoDriverNew(Domain *);
extern Rom::DriverInterface *elementSamplingDriverNew(Domain *);
extern Rom::DriverInterface *distrROMPostProcessingDriverNew(Domain *);
extern Rom::DriverInterface *distrElementSamplingDriverNew(Domain *);
extern Rom::DriverInterface *snapshotProjectionDriverNew(Domain *);
extern Rom::DriverInterface *positiveDualBasisDriverNew(Domain *);

#endif /* ROM_DRIVERINTERFACE_H */
