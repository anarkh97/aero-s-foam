#include <Rom.d/SnapshotNonLinDynamic.h>

#include <Corotational.d/GeomState.h>

SnapshotNonLinDynamic::SnapshotNonLinDynamic(Domain * domain) :
  NonLinDynamic(domain)
{}

void
SnapshotNonLinDynamic::saveStateSnapshot(const GeomState & snap) {
  // snap->extract(buffer);
  // TODO
}

void
SnapshotNonLinDynamic::saveResidualSnapshot(const Vector & snap) {
  // TODO
}
