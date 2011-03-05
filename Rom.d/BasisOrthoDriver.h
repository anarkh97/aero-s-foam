#ifndef ROM_BASIS_ORTHODRIVER_H
#define ROM_BASIS_ORTHODRIVER_H

#include "DriverInterface.h"

namespace Rom {

class BasisOrthoDriver : public DriverInterface {
public:
  virtual void solve();
  
  explicit BasisOrthoDriver(Domain *);

private:
  void preProcess();

  Domain *domain_;
};

} /* end namespace Rom */

Rom::DriverInterface *basisOrthoDriverNew(Domain *);

#endif /* ROM_BASIS_ORTHODRIVER_H */
