#ifndef PITA_HTS_NLLOCALNETWORK_H
#define PITA_HTS_NLLOCALNETWORK_H

#include "Fwk.h"
#include "Types.h"

#include "SliceMapping.h"

#include "../Seed.h"
#include "../RemoteState.h"

#include "HalfTimeSlice.h"
#include "CorrectionTimeSlice.h"

#include "../SeedErrorEvaluator.h"

#include <map>
#include <deque>

namespace Pita { namespace Hts {

class NlLocalNetwork : public Fwk::PtrInterface<NlLocalNetwork> {
public:
  EXPORT_PTRINTERFACE_TYPES(NlLocalNetwork);
private:
  DISALLOW_COPY_AND_ASSIGN(NlLocalNetwork);
};

} /* end namespace Hts */ } /* end namespace Pita */

#endif /* PITA_HTS_NLLOCALNETWORK_H */
