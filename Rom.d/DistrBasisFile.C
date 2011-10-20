#include "DistrBasisFile.h"

#include "DistrNodeDof6Buffer.h"

#include <Comm.d/Communicator.h>

#include <sys/types.h>

namespace Rom {

inline
MPI_Comm
DistrBasisOutputFile::nativeComm() const {
  return *comm_->getCommunicator();
}

inline
char *
DistrBasisOutputFile::getFormat() {
  return const_cast<char *>("native");
}

DistrBasisOutputFile::DistrBasisOutputFile(const std::string &fileName, int nodeCount, Communicator *comm) :
  fileName_(fileName),
  nodeCount_(nodeCount),
  stateCount_(0),
  status_(OUTDATED),
  comm_(comm)
{
  MPI_File_open(nativeComm(),
                const_cast<char *>(fileName_.c_str()),
                MPI_MODE_WRONLY | MPI_MODE_CREATE,
                MPI_INFO_NULL,
                &handle_);
  MPI_Offset currDisp = 0;
  MPI_File_set_size(handle_, currDisp);

  // Format
  MPI_File_set_view(handle_, currDisp, MPI_INT, MPI_INT, getFormat(), MPI_INFO_NULL);
  if (comm_->myID() == 0) {
    int swapFlag = 1;
    MPI_Status status;
    MPI_File_write(handle_, &swapFlag, 1, MPI_INT, &status);
  }

  MPI_File_get_size(handle_, &currDisp);
  MPI_File_set_view(handle_, currDisp, MPI_DOUBLE, MPI_DOUBLE, getFormat(), MPI_INFO_NULL); 
  
  if (comm_->myID() == 0) {
    double version = 1.0;
    MPI_Status status;
    MPI_File_write(handle_, &version, 1, MPI_DOUBLE, &status);
  }
  
  MPI_File_get_size(handle_, &currDisp);

  // Info
  const MPI_Datatype INTEGER_TYPE = MPI_UNSIGNED_LONG_LONG; 
  infoDisp_ = currDisp; 
  MPI_File_set_view(handle_, currDisp, INTEGER_TYPE, INTEGER_TYPE, getFormat(), MPI_INFO_NULL); 
  if (comm_->myID() == 0) {
    typedef u_int64_t BinIntType;
    BinIntType intBuffer[2];
    intBuffer[0] = static_cast<BinIntType>(stateCount_);
    intBuffer[1] = static_cast<BinIntType>(nodeCount_);
    
    MPI_Status status;
    MPI_File_write(handle_, &intBuffer, 2, INTEGER_TYPE, &status);
  }

  synchronizeFile(); 
  MPI_File_get_size(handle_, &dataDisp_);
  
  // Data
  MPI_File_set_view(handle_, dataDisp_, MPI_DOUBLE, MPI_DOUBLE, getFormat(), MPI_INFO_NULL); 
}

DistrBasisOutputFile::~DistrBasisOutputFile() {
  updateStateCountStatus(); // TODO: swallow exceptions
  MPI_File_close(&handle_);
}

void
DistrBasisOutputFile::stateAdd(const DistrNodeDof6Buffer &data) {
  const double defaultHeader = static_cast<double>(stateCount() + 1);
  stateAdd(data,defaultHeader);
}

void
DistrBasisOutputFile::stateAdd(const DistrNodeDof6Buffer &data, double headValue) {
  status_ = OUTDATED;

  MPI_Status status;
  
  // State header
  const MPI_Offset headerOffset = (nodeCount_ * 6 + 1) * stateCount_;
  if (comm_->myID() == 0) {
    MPI_File_write_at(handle_, headerOffset, &headValue, 1, MPI_DOUBLE, &status);
  }

  // State values 
  const MPI_Offset stateOffset = headerOffset + 1;
  for (DistrNodeDof6Buffer::NodeItConst it(data.globalNodeIndexBegin()),
                                        it_end(data.globalNodeIndexEnd());
       it != it_end; ++it) {
    const int iNode = *it;
    
    const MPI_Offset nodeOffset = stateOffset + (iNode * 6);
    MPI_File_write_at(handle_, nodeOffset, const_cast<double *>(data[iNode]), 6, MPI_DOUBLE, &status);
  }
  
  ++stateCount_;
}

void
DistrBasisOutputFile::updateStateCountStatus() {
  synchronizeFile();

  const MPI_Datatype INTEGER_TYPE = MPI_UNSIGNED_LONG_LONG; 
  MPI_File_set_view(handle_, infoDisp_, INTEGER_TYPE, INTEGER_TYPE, getFormat(), MPI_INFO_NULL); 
  if (comm_->myID() == 0) {
    typedef u_int64_t BinIntType;
    BinIntType intBuffer;
    intBuffer = static_cast<BinIntType>(stateCount_);
    
    MPI_Status status;
    MPI_File_write_at(handle_, 0, &intBuffer, 1, INTEGER_TYPE, &status);
  }
  
  synchronizeFile();
  status_ = UP_TO_DATE;
  MPI_File_set_view(handle_, dataDisp_, MPI_DOUBLE, MPI_DOUBLE, getFormat(), MPI_INFO_NULL); 
}

DistrBasisInputFile::DistrBasisInputFile(const std::string &fileName) :
  BasisBinaryInputFile(fileName)
{}

void
DistrBasisOutputFile::synchronizeFile() {
  MPI_File_sync(handle_);
  MPI_Barrier(nativeComm());
  MPI_File_sync(handle_);
}

const DistrNodeDof6Buffer&
DistrBasisInputFile::currentStateBuffer(DistrNodeDof6Buffer &target) {
  if (validCurrentState()) {
    seekCurrentState();

    typedef DistrNodeDof6Buffer::NodeItConst NodeIt;
    const NodeIt it_end = target.globalNodeIndexEnd();
    for (NodeIt it = target.globalNodeIndexBegin(); it != it_end; ++it) {
      const int iNode = *it;
      seekNode(iNode);
      readCurrentNode(target[iNode]);
    }

    seekNode(nodeCount());
    validateNextOffset();
  }

  return target;
}


} /* end namespace Rom */
