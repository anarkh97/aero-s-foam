#ifndef ROM_DISTRELEMENTSAMPLINGDRIVER_H
#define ROM_DISTRELEMENTSAMPLINGDRIVER_H

#include "DriverInterface.h"
#include <memory>
#include "BasisId.h"

class Communicator;
class DistrInfo;
class Connectivity;

namespace Rom {

template <typename Scalar, template <typename> class GenVecType> class GenVecBasis;
typedef GenVecBasis<double, GenDistrVector> DistrVecBasis;

class DistrElementSamplingDriver :public MultiDomainDynam, public DriverInterface {
public:
  virtual void solve();
  
  DistrElementSamplingDriver(Domain *, Communicator *);
  const DistrInfo& vectorSize() const;

private:
  Communicator *comm_;
  std::auto_ptr<Connectivity> elemToNode_;
  int snapSize(BasisId::Type type);
  void buildDomainCdsa();
};

} /* end namespace Rom */

Rom::DriverInterface *distrElementSamplingDriverNew(Domain *);

#endif /* ROM_DISTRELEMENTSAMPLINGDRIVER_H */
