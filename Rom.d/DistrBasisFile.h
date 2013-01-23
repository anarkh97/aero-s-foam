#ifndef ROM_DISTRBASISFILE_H
#define ROM_DISTRBASISFILE_H

#include "BasisBinaryFile.h"

#include "NodeDof6Buffer.h"
#include "DistrNodeDof6Buffer.h"

#include <Comm.d/Communicator.h>

#include <string>
#include <vector>
#include <map>
#include <iterator>
#include <numeric>
#include <memory>
#include <cstddef>

#include <cassert>

namespace Rom {

class DistrBasisOutputFile : public BasisBinaryFile {
public:
  const std::string &fileName() const { return binFile_->pathName(); }

  int nodeCount() const  { return binFile_->itemCount(); }
  int stateCount() const { return binFile_->stateCount(); }

  void stateAdd(const DistrNodeDof6Buffer &data);
  void stateAdd(const DistrNodeDof6Buffer &data, double headValue);

  // Constructor
  // Collective operation (must be called by all processes sharing the communicator)
  template <typename IdxIt>
  DistrBasisOutputFile(const std::string &fileName, int globalNodeCount, IdxIt localIdxBegin, IdxIt localIdxEnd, Communicator *comm);

private:
  template <typename IdxIt>
  void
  resetHandler(const std::string &fileName, int globalNodeCount, int localOffset, IdxIt localIdxBegin, IdxIt localIdxEnd) {
    binFile_.reset(new BinaryResultOutputFile(fileName, NODAL_DATA_FLAG, DESC, globalNodeCount, DOFS_PER_NODE, localOffset, localIdxBegin, localIdxEnd, VERSION));
  }

  std::auto_ptr<BinaryResultOutputFile> binFile_;

  // Disallow copy and assigment
  DistrBasisOutputFile(const DistrBasisOutputFile &);
  DistrBasisOutputFile &operator=(const DistrBasisOutputFile &);
};

template <typename Scalar>
Scalar
distr_exclusive_partial_sum(Scalar v, Communicator *comm) {
  const int cpuCount = comm->numCPUs();
  std::vector<Scalar> a(cpuCount, Scalar());
  
  typename std::vector<Scalar>::iterator it(a.begin() + comm->myID());
  *it = v;
  comm->globalSum(a.size(), &a[0]);

  return std::accumulate(a.begin(), it, Scalar());
}

template <typename Scalar, typename Scalar2>
Scalar
asserting_distr_exclusive_partial_sum(Scalar v, Communicator *comm, Scalar2 sum) {
  const int cpuCount = comm->numCPUs();
  std::vector<Scalar> a(cpuCount, Scalar());
  
  typename std::vector<Scalar>::iterator it(a.begin() + comm->myID());
  *it = v;
  comm->globalSum(a.size(), &a[0]);

  assert(std::accumulate(a.begin(), a.end(), Scalar()) == sum);

  return std::accumulate(a.begin(), it, Scalar());
}

template <typename IdxIt>
DistrBasisOutputFile::DistrBasisOutputFile(const std::string &fileName, int globalNodeCount, IdxIt localIdxBegin, IdxIt localIdxEnd, Communicator *comm) :
  binFile_(NULL)
{
  const int localOffset = asserting_distr_exclusive_partial_sum(std::distance(localIdxBegin, localIdxEnd), comm, globalNodeCount);
  if (localOffset == 0) {
    // Master process goes first, creates/truncates the file without interference
    resetHandler(fileName, globalNodeCount, localOffset, localIdxBegin, localIdxEnd);
    comm->sync();
  } else {
    // Other processes wait until master process has created the file
    comm->sync();
    resetHandler(fileName, globalNodeCount, localOffset, localIdxBegin, localIdxEnd);
  }
#if defined(_POSIX_THREAD_SAFE_FUNCTIONS) && defined(_AEROS_ASYCHRONOUS_IO)
  // When this is macro defined we can safely assume the fwrite is thread-safe
  // however, since we are using fopen and there is no way to open a file for use with fseek/fwrite without truncating
  // we need another sync to guarantee thread safety
  // see comments in Utils.d/BinaryResultFile.h
  comm->sync();
  // now the Master process can call writePrelude
  if (localOffset == 0) {
    binFile_->writePrelude();
  }
#endif
}

class DistrBasisInputFile : public BasisBinaryInputFile {
public:
  explicit DistrBasisInputFile(const std::string &fileName);
  
  const DistrNodeDof6Buffer &currentStateBuffer(DistrNodeDof6Buffer &target);
  using BasisBinaryInputFile::currentStateBuffer; // Do not hide inherited member function

private:
  std::map<int, int> fileNodeIds_;
  NodeDof6Buffer fileBuffer_; 
};

} /* end namespace Rom */

#endif /* ROM_DISTRBASISFILE_H */
