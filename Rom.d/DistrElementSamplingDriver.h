#ifndef ROM_DISTRELEMENTSAMPLINGDRIVER_H
#define ROM_DISTRELEMENTSAMPLINGDRIVER_H

#include "DriverInterface.h"
#include <memory>
#include "BasisId.h"

class Communicator;
class DistrInfo;

namespace Rom {

class DistrElementSamplingDriver : public MultiDomainDynam, public DriverInterface {
public:
  virtual void solve();
  
  DistrElementSamplingDriver(Domain *, Communicator *);
  const DistrInfo& vectorSize() const;

private:
  Communicator *comm_;
  void buildDomainCdsa();
  void subMakeMass(int isub, SparseMatrix **subM);
};

} /* end namespace Rom */

Rom::DriverInterface *distrElementSamplingDriverNew(Domain *);

#endif /* ROM_DISTRELEMENTSAMPLINGDRIVER_H */
