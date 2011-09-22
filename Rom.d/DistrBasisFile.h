#ifndef ROM_DISTRBASISFILE_H
#define ROM_DISTRBASISFILE_H

#include "BasisBinaryFile.h"

#include <mpi.h>

class Communicator;

#include <string>

namespace Rom {

class DistrNodeDof6Buffer;

class DistrBasisOutputFile {
public:
  const std::string &fileName() const { return fileName_; }

  int nodeCount() const  { return nodeCount_;  }
  int stateCount() const { return stateCount_; }

  enum StateCountStatus { UP_TO_DATE, OUTDATED };
  StateCountStatus stateCountStatus() const { return status_; }
  void updateStateCountStatus();

  void stateAdd(const DistrNodeDof6Buffer &data);
  void stateAdd(const DistrNodeDof6Buffer &data, double headValue);

  DistrBasisOutputFile(const std::string &fileName, int nodeCount, Communicator *comm);

  ~DistrBasisOutputFile();

private:
  MPI_Comm nativeComm() const;

  void synchronizeFile();

  const std::string fileName_;
  const int nodeCount_;
  int stateCount_;

  StateCountStatus status_;

  Communicator *comm_;
  MPI_File handle_;
  MPI_Offset infoDisp_;
  MPI_Offset dataDisp_;
};

class DistrBasisInputFile : public BasisBinaryInputFile {
public:
  explicit DistrBasisInputFile(const std::string &fileName);

  // Hides BasisBinaryInputFile::currentStateBuffer
  const DistrNodeDof6Buffer &currentStateBuffer(DistrNodeDof6Buffer &target);
};

} /* end namespace Rom */

#endif /* ROM_DISTRBASISFILE_H */
