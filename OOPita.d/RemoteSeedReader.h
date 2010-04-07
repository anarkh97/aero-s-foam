#ifndef PITA_REMOTESEEDREADER_H
#define PITA_REMOTESEEDREADER_H

#include "Fwk.h"
#include "Types.h"

#include "Seed.h"

#include "RemoteSharedState.h"

namespace Pita {

typedef RemoteSharedStateReader<DynamState> RemoteSeedReader;

template <>
void RemoteSharedStateReader<DynamState>::statusIs(Status s);

} // end namespace Pita

#endif /* PITA_REMOTESEEDREADER_H */
