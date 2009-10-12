#ifndef PITA_COMMSPLITTER_H
#define PITA_COMMSPLITTER_H

#include "Fwk.h"
#include "Types.h"

#include <Comm.d/Communicator.h>

namespace Pita {

class CommSplitter : public Fwk::PtrInterface<CommSplitter> {
public:
  EXPORT_PTRINTERFACE_TYPES(CommSplitter);

  enum CpuColor {
    FINE   = 0,
    COARSE = 1
  };

  CpuColor localColor() const { return localColor_; }
  CpuRank coarseCpu() const { return defaultCoarseCpu_; }

  Communicator * originComm() const { return originComm_; } // Not owned by CommSplitter
  Communicator * splitComm() const { return splitComm_; } // Owned by CommSplitter
  Communicator * interComm() const { return interComm_; } // Owned by CommSplitter

  static Ptr New(Communicator * originComm) {
    return new CommSplitter(originComm);
  }

protected:
  explicit CommSplitter(Communicator * originComm); // Must remain valid while CommSplitter is alive

  ~CommSplitter();

  static const CpuRank defaultCoarseCpu_;
  static const CpuRank defaultFineLeader_;

private:
  CpuRank coarseCpu_;
  CpuColor localColor_;
  Communicator * originComm_;
  Communicator * splitComm_;
  Communicator * interComm_;
};

} /* end namespace Pita */

#endif /* PITA_COMMSPLITTER_H */
