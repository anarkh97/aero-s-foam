#ifndef PITA_HTS_NLBASISUPDATE_H
#define PITA_HTS_NLBASISUPDATE_H

#include "Fwk.h"
#include "Types.h"

#include "GlobalStateSharing.h"

namespace Pita { namespace Hts {

class NlBasisUpdate : public Fwk::PtrInterface<NlBasisUpdate> {
public:
  EXPORT_PTRINTERFACE_TYPES(NlBasisUpdate);

  GlobalStateSharing * globalUpdate() { return globalUpdate_.ptr(); }

  static Ptr New(Communicator * timeComm, size_t vectorSize, GlobalStateSharing::Strategy globalStrategy) {
    return new NlBasisUpdate(timeComm, vectorSize, globalStrategy);
  }

protected:
  NlBasisUpdate(Communicator *, size_t, GlobalStateSharing::Strategy);

private:
  GlobalStateSharing::Ptr globalUpdate_;

  DISALLOW_COPY_AND_ASSIGN(NlBasisUpdate);
};

} /* end namespace Hts */ } /* end namespace Pita */

#endif /* PITA_HTS_NLBASISUPDATE_H */
