#include "BasisBinaryFile.h"

#include <sys/types.h>

#include <stdexcept>
#include <cassert>
#include <iostream>

namespace Rom {

const double BasisBinaryFile::VERSION = 2.0;
const std::string BasisBinaryFile::DESC = "rob";
const int BasisBinaryFile::NODAL_DATA_FLAG = 1;
const int BasisBinaryFile::DOFS_PER_NODE = 6;


BasisBinaryOutputFile::BasisBinaryOutputFile(const std::string &fileName, int nodeCount, bool restart) :
  binFile_(fileName, NODAL_DATA_FLAG, DESC, nodeCount, DOFS_PER_NODE, VERSION, restart)
{}

void
BasisBinaryOutputFile::stateAdd(const NodeDof6Buffer &data, double headValue) {
  assert(nodeCount() == data.size());

  // Dump all information in one pass
  binFile_.stateAdd(headValue, data.array());
}


BasisBinaryInputFile::BasisBinaryInputFile(const std::string &fileName) :
  binFile_(fileName)
{

//  fprintf(stderr," binFile_.version = %f matches VERSION %f for file ", binFile_.version(),VERSION);
//  std::cerr << fileName << std::endl;

//  fprintf(stderr," binFile_.dataType = %d and NODAL_DATA_FLAG = %d \n", binFile_.dataType(), NODAL_DATA_FLAG);
//  fprintf(stderr," binFile_.description = %s and DESC = %s \n", binFile_.description(), DESC);
//  fprintf(stderr," binFile_.itemDimension = %d and DOFS_PER_NODE = %d \n", binFile_.itemDimension(), DOFS_PER_NODE);

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
  assert(nodeCount() == target.size());
  
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
