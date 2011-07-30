#include "BasisBinaryFile.h"

#include <sys/types.h>

#include <stdexcept>

namespace Rom {

namespace {
  const double VERSION = 1.0;
  const char *WRITE_FLAG = "w";
  const char *READ_FLAG = "r";
  typedef u_int64_t BinIntType;
  typedef BinFileHandler::OffType OffSetType;
} // end anonymous namespace

BasisBinaryOutputFile::BasisBinaryOutputFile(const std::string &fileName, int nodeCount) :
  fileName_(fileName),
  nodeCount_(nodeCount),
  stateCount_(0),
  binHandler_(fileName.c_str(), WRITE_FLAG, VERSION) 
{
  // WARNING: The program exits when it fails to open the file
  // TODO: throw std::runtime_error instead

  BinIntType intBuffer[2];
  intBuffer[0] = static_cast<BinIntType>(stateCount_);
  intBuffer[1] = static_cast<BinIntType>(nodeCount_);
  binHandler_.write(intBuffer, 2);
}

BasisBinaryOutputFile::StateCountStatus
BasisBinaryOutputFile::stateCountStatus() const {
  return UP_TO_DATE;
}

void
BasisBinaryOutputFile::updateStateCountStatus() {
  // Nothing to do: State count automatically kept up-to-date
}

void
BasisBinaryOutputFile::writeStateHeader(double headValue) {
  binHandler_.write(&headValue, 1);
}

void
BasisBinaryOutputFile::incrementStateCount() {
  const int newStateCount = stateCount_ + 1;

  // Save the current offset before rewinding
  const OffSetType offset = binHandler_.tell(); 
  binHandler_.seek(OffSetType(0));

  // Write the state count in a binary compatible format
  const BinIntType binStateCount = static_cast<BinIntType>(newStateCount);
  binHandler_.write(&binStateCount, 1);

  // Restore the offset
  binHandler_.seek(offset);

  // Commit changes
  stateCount_ = newStateCount;
}

BasisBinaryInputFile::BasisBinaryInputFile(const std::string &fileName) :
  fileName_(fileName),
  nodeCount_(),
  stateCount_(),
  currentStateIndex_(),
  currentStateHeaderValue_(),
  binHandler_(fileName.c_str(), READ_FLAG, VERSION),
  savedOffset_(),
  validNextOffset_(false)
{
  // WARNING: The program exits when it fails to open the file
  // TODO: throw std::runtime_error instead
  
  if (binHandler_.getVersion() != VERSION) {
    throw std::runtime_error("Incompatible binary file version");
  }

  // Read file header
  BinIntType intBuffer[2];
  binHandler_.read(intBuffer, 2);
  stateCount_ = static_cast<int>(intBuffer[0]);
  nodeCount_ = static_cast<int>(intBuffer[1]);
  validNextOffset_ = true;

  if (validCurrentState()) {
    seekNextState();
  }
}

void
BasisBinaryInputFile::currentStateIndexInc() {
  if (validCurrentState()) {
    seekNextState();
    ++currentStateIndex_; 
  }
}

void
BasisBinaryInputFile::seekNextState() {
  if (!validNextOffset_) {
    throw std::runtime_error("Must read current state before advancing");
  }
  binHandler_.read(&currentStateHeaderValue_, 1);
  savedOffset_ = binHandler_.tell();
  validNextOffset_ = false;
}

void
BasisBinaryInputFile::seekCurrentState() {
  if (validNextOffset_) {
    binHandler_.seek(savedOffset_);
  }
}

} /* end namespace Rom */
