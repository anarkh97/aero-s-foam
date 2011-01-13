#include "BasisOutputFile.h"

#include <string>
#include <cstdio>

#include <stdexcept>

#include <cassert>

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
BasisOutputFile::stateAdd(const double(*data)[6]) {
  stateAdd(data, static_cast<double>(stateCount_ + 1));
}

void
BasisOutputFile::stateAdd(const double(*data)[6], double headValue) {
  const int w = width_;
  const int p = precision_;

  std::fprintf(stream_, "  % *.*E\n", w, p, headValue);
  
  for (int iNode = 0; iNode < nodeCount_; iNode++)  {
    std::fprintf(stream_, " % *.*E % *.*E % *.*E % *.*E % *.*E % *.*E\n",
        w, p, data[iNode][0], w, p, data[iNode][1], w, p, data[iNode][2],
        w, p, data[iNode][3], w, p, data[iNode][4], w, p, data[iNode][5]);
  }

  stateCount_++;
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
