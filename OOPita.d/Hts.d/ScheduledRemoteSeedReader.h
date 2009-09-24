#ifndef PITA_HTS_SCHEDULEDREMOTESEEDREADER_H
#define PITA_HTS_SCHEDULEDREMOTESEEDREADER_H

#include "../RemoteSeedReader.h"
#include "ScheduledRemoteSharedState.h"

namespace Pita {

typedef ScheduledRemoteSharedStateReader<DynamState> ScheduledRemoteSeedReader;

} // end namespace Pita

#endif /* PITA_HTS_SCHEDULEDREMOTESEEDREADER_H */