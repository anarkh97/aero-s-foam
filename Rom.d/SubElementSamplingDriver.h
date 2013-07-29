#ifndef ROM_SUBELEMENTSAMPLINGDRIVER_H
#define ROM_SUBELEMENTSAMPLINGDRIVER_H

#include "ElementSamplingDriver.h"

namespace Rom {

class SubElementSamplingDriver : public ElementSamplingDriver<std::vector<double>,size_t> {
public:
  explicit SubElementSamplingDriver(Domain *);
  void getGlobalWeights(Vector &solution, vector<double> &lweights, vector<int> &lelemids, bool verboseFlag = true);
private:
  void preProcess();
};

} // end namespace Rom

Rom::DriverInterface *subElementSamplingDriverNew(Domain *);

#endif /* ROM_SUBELEMENTSAMPLINGDRIVER_H */
