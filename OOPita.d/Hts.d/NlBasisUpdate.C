#include "NlBasisUpdate.h"

namespace Pita { namespace Hts {

NlBasisUpdate::NlBasisUpdate(Communicator * timeComm, size_t vectorSize, GlobalStateSharing::Strategy globalStrategy) :
  globalUpdate_(new GlobalStateSharing(timeComm, vectorSize, globalStrategy))
{}

} /* end namespace Hts */ } /* end namespace Pita */
