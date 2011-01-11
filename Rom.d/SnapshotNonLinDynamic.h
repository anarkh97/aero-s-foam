#ifndef ROM_SNAPSHOTNONLINDYNAMIC_H
#define ROM_SNAPSHOTNONLINDYNAMIC_H

#include <Problems.d/NonLinDynam.h>

class Domain;
class GeomState;

class SnapshotNonLinDynamic : public NonLinDynamic {
public:
  SnapshotNonLinDynamic(Domain *);

  void saveStateSnapshot(const GeomState &);
  void saveResidualSnapshot(const Vector &);
};

#include "Fom.h"

#endif /* ROM_SNAPSHOTNONLINDYNAMIC_H */
