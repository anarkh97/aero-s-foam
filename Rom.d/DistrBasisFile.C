#include "DistrBasisFile.h"

#include "RenumberingUtils.h"

#include <algorithm>
#include <iterator>
#include <stdexcept>

#include <cassert>

namespace Rom {

void
DistrBasisOutputFile::stateAdd(const DistrNodeDof6Buffer &data) {
  const double defaultHeader = static_cast<double>(stateCount() + 1);
  stateAdd(data,defaultHeader);
}

void
DistrBasisOutputFile::stateAdd(const DistrNodeDof6Buffer &data, double headValue) {
  assert(data.localNodeCount() == binFile_->localItemCount());
  SimpleBuffer<double> buffer(binFile_->localDataSize());

  double *target = buffer.array();
  for (BinaryResultOutputFile::ItemIdIterator it = binFile_->itemIdBegin(),
                                              it_end = binFile_->itemIdEnd();
                                              it != it_end;
                                              ++it) {
    const double *origin = data[*it];
    std::copy(origin, origin + DOFS_PER_NODE, target);
    target += DOFS_PER_NODE;
  }

  binFile_->stateAdd(headValue, buffer.array());
}


DistrBasisInputFile::DistrBasisInputFile(const std::string &fileName) :
  BasisBinaryInputFile(fileName),
  fileNodeIds_(),
  fileBuffer_(nodeCount())
{
  inverse_numbering(nodeIdBegin(), nodeIdEnd(), std::inserter(fileNodeIds_, fileNodeIds_.end()));
  assert(fileNodeIds_.size() == nodeCount());
}

const DistrNodeDof6Buffer &
DistrBasisInputFile::currentStateBuffer(DistrNodeDof6Buffer &target) {
  // Retrieve all information in the internal buffer
  // TODO: More economical approach
  BasisBinaryInputFile::currentStateBuffer(fileBuffer_);
  
  typedef DistrNodeDof6Buffer::NodeItConst NodeIt;
  const NodeIt it_end = target.globalNodeIndexEnd();
  for (NodeIt it = target.globalNodeIndexBegin(); it != it_end; ++it) {
    const int iNode = *it; // Id of the node requested by the local buffer

    // Location of requested node in internal buffer
    std::map<int, int>::const_iterator it_loc = fileNodeIds_.find(iNode);
    if (it_loc == fileNodeIds_.end()) {
      throw std::runtime_error("Requested nodal data missing from file");
    }
    const double *nodeBuffer = fileBuffer_[it_loc->second];

    std::copy(nodeBuffer, nodeBuffer + DOFS_PER_NODE, &target[iNode][0]);
  }

  return target;
}

} /* end namespace Rom */
