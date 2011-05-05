#include "BasisOutputFile.h"

#include <string>
#include <cstdio>

#include <stdexcept>

#include <cassert>

namespace Rom {

const int
BasisOutputFile::STATE_COUNT_LENGTH = 10;

BasisOutputFile::BasisOutputFile(const std::string &fileName, int nodeCount) :
  fileName_(fileName),
  nodeCount_(nodeCount),
  width_(23),
  precision_(15),
  stateCount_(0),
  stream_(NULL),
  stateCountOnFile_(0)
{
  stream_ = std::fopen(fileName_.c_str(), "wt");

  if (!stream_) {
   throw std::runtime_error("Cannot open output file"); 
  }

  writeStateCount();
  writeNodeCount();
}

BasisOutputFile::~BasisOutputFile() {
  if (stateCountStatus() != UP_TO_DATE) {
    try {
      rewindAndWriteStateCount();
    } catch (std::runtime_error &) {
      // Log error and swallow exception
      fprintf(stderr, "WARNING: Corrupted output file %s\n", fileName_.c_str());
    }
  }

  std::fclose(stream_);
}

void
BasisOutputFile::updateStateCountStatus() {
  if (stateCountStatus() == UP_TO_DATE) {
    return;
  }

  rewindAndWriteStateCount();

  std::fflush(stream_);
  std::fseek(stream_, 0, SEEK_END);
}

void
BasisOutputFile::rewindAndWriteStateCount() {
  std::fflush(stream_);
  std::rewind(stream_);

  writeStateCount();
}

void
BasisOutputFile::writeStateCount() {
  const int actualLength = std::fprintf(stream_, "%-*d", STATE_COUNT_LENGTH, stateCount_);
  if (actualLength != STATE_COUNT_LENGTH) {
    throw std::runtime_error("Corrupted state count");
  }
  std::fprintf(stream_, "\n");
  stateCountOnFile_ = stateCount_;
}

void
BasisOutputFile::writeNodeCount() {
  fprintf(stream_, "%d\n", nodeCount_);
}

void
BasisOutputFile::writeStateHeader(double value) {
  std::fprintf(stream_, "  % *.*E\n", width_, precision_, value);
}

} /* end namespace Rom */
