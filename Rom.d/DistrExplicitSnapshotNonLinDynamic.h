#ifndef ROM_DISTREXPLICITSNAPSHOTNONLINDYNAMIC_H
#define ROM_DISTREXPLICITSNAPSHOTNONLINDYNAMIC_H

#include <Paral.d/MDDynam.h>

#include <Problems.d/DynamProbTraits.h>

namespace Rom {

class DistrExplicitSnapshotNonLinDynamic : public MultiDomainDynam {
public:
  explicit DistrExplicitSnapshotNonLinDynamic(Domain *);

  // Overriding via hiding
  void preProcess(); // Required additional pre-processing

  // Added functionality
  void snapshotAdd(const DistrVector &s);

  ~DistrExplicitSnapshotNonLinDynamic();

private:
  class SnapshotHandler;
  friend class SnapshotHandler;
  SnapshotHandler * snapshotHandler_;

  // Disallow copy and assignment
  DistrExplicitSnapshotNonLinDynamic(const DistrExplicitSnapshotNonLinDynamic &);
  DistrExplicitSnapshotNonLinDynamic &operator=(const DistrExplicitSnapshotNonLinDynamic &);
};

} // end namespace Rom

inline
void
handleDisplacement(Rom::DistrExplicitSnapshotNonLinDynamic &probDesc, DistrVector &d) {
  probDesc.snapshotAdd(d);
}

#endif /* ROM_DISTREXPLICITSNAPSHOTNONLINDYNAMIC_H */
