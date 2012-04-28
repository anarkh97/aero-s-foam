#ifndef ROM_DISTREXPLICITPODPROJECTIONNONLINDYNAMIC_H
#define ROM_DISTREXPLICITPODPROJECTIONNONLINDYNAMIC_H

#include "DistrExplicitPodProjectionNonLinDynamicBase.h"

#include <Problems.d/DynamProbTraits.h>

#include <memory>

namespace Rom {

class DistrExplicitPodProjectionNonLinDynamic : public DistrExplicitPodProjectionNonLinDynamicBase {
public:
  explicit DistrExplicitPodProjectionNonLinDynamic(Domain *);

  // Overriding via hiding
  void preProcess(); // Additional pre-processing

  // Added functionality
  void forceSnapshotAdd(const DistrVector &s);

  ~DistrExplicitPodProjectionNonLinDynamic();

protected:
  class SnapshotHandler;

private:
  friend class SnapshotHandler;
  std::auto_ptr<SnapshotHandler> snapshotHandler_;
};

} // end namespace Rom

inline
void
handleForce(Rom::DistrExplicitPodProjectionNonLinDynamic &probDesc, DistrVector &f) {
  probDesc.forceSnapshotAdd(f);
}

#endif /* ROM_DISTREXPLICITPODPROJECTIONNONLINDYNAMIC_H */
