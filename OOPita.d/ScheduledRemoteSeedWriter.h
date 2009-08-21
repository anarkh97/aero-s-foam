#ifndef PITA_SCHEDULEDREMOTESEEDWRITER_H
#define PITA_SCHEDULEDREMOTESEEDWRITER_H

#include "RemoteSeedWriter.h"

#include "ScheduledRemoteSharedState.h"
#include "DynamState.h"

namespace Pita {

typedef ScheduledRemoteSharedStateWriter<DynamState> ScheduledRemoteSeedWriter;

} // end namespace Pita

#endif /* PITA_SCHEDULEDREMOTESEEDWRITER_H */
