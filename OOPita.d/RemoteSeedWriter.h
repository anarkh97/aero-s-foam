#ifndef PITA_REMOTESEEDWRITER_H
#define PITA_REMOTESEEDWRITER_H

#include "Fwk.h"
#include "Types.h"
#include "HalfSliceTypes.h"

#include "RemoteSharedState.h"
#include "Seed.h"

namespace Pita {

typedef RemoteSharedStateWriter<DynamState> RemoteSeedWriter;

template<>
void RemoteSharedStateWriter<DynamState>::statusIs(Status s);

} // end namespace Pita

#endif
