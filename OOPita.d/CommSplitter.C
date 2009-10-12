#include "CommSplitter.h"

namespace Pita {

const CpuRank CommSplitter::defaultCoarseCpu_ = CpuRank(0);
const CpuRank CommSplitter::defaultFineLeader_ = CpuRank(1);

CommSplitter::CommSplitter(Communicator * originComm) :
  localColor_(FINE),
  originComm_(originComm),
  splitComm_(NULL),
  interComm_(NULL)
{
  // Determine local task
  int originId = originComm->myID();
  if (CpuRank(originComm->myID()) == defaultCoarseCpu_) {
    localColor_ = COARSE;
  }

  // Create non-overlapping communicators
  int originCpus = originComm->numCPUs();
  MPI_Comm * originMComm = originComm->getCommunicator();

  MPI_Comm splitMComm;
  MPI_Comm_split(*originMComm, localColor_, originId, &splitMComm);

  // Create intercommunicator
  int splitCpus, splitId;
  MPI_Comm_size(splitMComm, &splitCpus);
  MPI_Comm_rank(splitMComm, &splitId); 

  MPI_Comm interMComm;
  CpuRank remoteLeader = (localColor_ == COARSE) ? defaultFineLeader_ : defaultCoarseCpu_;
  MPI_Intercomm_create(splitMComm, 0, *originMComm, remoteLeader.value(), 0, &interMComm);
 
  // Publish results
  splitComm_ = new Communicator(splitMComm, stderr);
  interComm_ = new Communicator(interMComm, stderr);
}

CommSplitter::~CommSplitter() {
  // Mark MPI Communicator for deletion
  MPI_Comm_free(interComm_->getCommunicator());
  MPI_Comm_free(splitComm_->getCommunicator());
  
  delete interComm_;
  delete splitComm_;
}

} /* end namespace Pita */
