#include "VecBasisFile.h"

#include "BasisFileStream.h"

#include <Math.d/Vector.h>

#include <algorithm>

BasisInputStream &
operator>>(BasisInputStream &in, VecBasis &sink) {
  GenVecBasis<double> temp(in.size() - in.currentVectorRank(), in.vectorSize());
  readVectors(in, temp.begin(), temp.end());
  sink.swap(temp);

  return in;
}

BasisInputStream &
readVectors(BasisInputStream &in, VecBasis &sink, int countMax) {
  const int count = std::max(std::min(in.size() - in.currentVectorRank(), countMax), 0);
  GenVecBasis<double> temp(count, in.vectorSize());
  readVectors(in, temp.begin(), temp.end());
  sink.swap(temp);

  return in;
}

BasisOutputStream &
operator<<(BasisOutputStream &out, const VecBasis &source) {
  return writeVectors(out, source.begin(), source.end());
}

BasisOutputStream &
writeVectors(BasisOutputStream &out, const VecBasis &source, int countMax) {
  VecBasis::const_iterator last = std::min(source.end(), source.begin() + countMax);
  return writeVectors(out, source.begin(), last);
}
