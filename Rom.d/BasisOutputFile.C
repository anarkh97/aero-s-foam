#include "BasisOutputFile.h"

#include <string>
#include <cstdio>

#include <stdexcept>

BasisOutputFile::BasisOutputFile(const std::string &fileName, int nodeCount):
  fileName_(fileName),
  nodeCount_(nodeCount),
  width_(23),
  precision_(15),
  stream_(NULL),
  stateCount_(0)
{
  stream_ = std::fopen(fileName_.c_str(), "wt");

  if (!stream_) {
   throw std::runtime_error("Cannot open output file"); 
  }
}

BasisOutputFile::~BasisOutputFile() {
  std::fclose(stream_);
}

void
BasisOutputFile::stateAdd(const double(*data)[6]) {
  const int w = width_;
  const int p = precision_;

  std::fprintf(stream_, "  % *.*E\n", w, p, static_cast<double>(stateCount_ + 1));
  
  for (int iNode = 0; iNode < nodeCount_; iNode++)  {
    std::fprintf(stream_, " % *.*E % *.*E % *.*E % *.*E % *.*E % *.*E\n",
        w, p, data[iNode][0], w, p, data[iNode][1], w, p, data[iNode][2],
        w, p, data[iNode][3], w, p, data[iNode][4], w, p, data[iNode][5]);
  }

  stateCount_++;
}
