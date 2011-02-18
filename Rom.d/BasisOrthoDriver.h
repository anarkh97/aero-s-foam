#ifndef ROM_BASIS_ORTHODRIVER_H
#define ROM_BASIS_ORTHODRIVER_H

#include "DriverInterface.h"

class BasisOrthoDriver : public RomDriverInterface {
public:
  virtual void solve();
  
  explicit BasisOrthoDriver(Domain *);

private:
  void preProcess();

  Domain *domain_;
};

RomDriverInterface *basisOrthoDriverNew(Domain *);

#endif /* ROM_BASIS_ORTHODRIVER_H */
