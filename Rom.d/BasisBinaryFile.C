#include "BasisBinaryFile.h"

#include <sys/types.h>

#include <stdexcept>
#include <cassert>

namespace Rom {

const double BasisBinaryFile::VERSION = 2.0;
const std::string BasisBinaryFile::DESC = "rob";
const int BasisBinaryFile::NODAL_DATA_FLAG = 1;
const int BasisBinaryFile::DOFS_PER_NODE = 6;


BasisBinaryOutputFile::BasisBinaryOutputFile(const std::string &fileName, int nodeCount) :
  binFile_(fileName, NODAL_DATA_FLAG, DESC, nodeCount, DOFS_PER_NODE, VERSION)
{}


BasisBinaryInputFile::BasisBinaryInputFile(const std::string &fileName) :
  binFile_(fileName)
{
  if (binFile_.version() != VERSION) {
    throw std::runtime_error("Incompatible binary file version");
  }

  if (binFile_.dataType() != NODAL_DATA_FLAG) {
    throw std::runtime_error("Non-nodal data");
  }
  
  if (binFile_.description() != DESC) {
    throw std::runtime_error("Incorrect description");
  }
  
  if (binFile_.itemDimension() != DOFS_PER_NODE) {
    throw std::runtime_error("Incorrect #dofs/node");
  }
  
  assert(currentStateIndex() == 0);

  cacheStateHeaderValue();
}

const NodeDof6Buffer &
BasisBinaryInputFile::currentStateBuffer(NodeDof6Buffer &target) {
  assert(validCurrentState());
  
  // Retrieve all information in one pass
  binFile_.state(target.array());

  return target;
}

void
BasisBinaryInputFile::currentStateIndexInc() {
  if (validCurrentState()) {
    binFile_.stateRankInc();
    cacheStateHeaderValue();
  }
}

void BasisBinaryInputFile::cacheStateHeaderValue() {
  if (validCurrentState()) {
    currentStateHeaderValue_ = binFile_.stateStamp();
  }
}

} /* end namespace Rom */
