#include "VecBasisFile.h"

#include "BasisFileStream.h"
#include "NodalRestrictionMapping.h"

#include <Math.d/Vector.h>

#include <cassert>
#include <algorithm>

namespace Rom {

// Output
BasisOutputStream &
operator<<(BasisOutputStream &out, const VecBasis &source) {
  assert(out.vectorSize() == source.vectorSize());
  return writeVectors(out, source.begin(), source.end());
}

BasisOutputStream &
writeVectors(BasisOutputStream &out, const VecBasis &source, int countMax) {
  assert(out.vectorSize() == source.vectorSize());
  VecBasis::const_iterator last = std::min(source.end(), source.begin() + countMax);
  return writeVectors(out, source.begin(), last);
}

BasisOutputStream &
writeExtendedVectors(BasisOutputStream &out, const VecBasis &source, const NodalRestrictionMapping &mapping) {
  assert(out.vectorSize() == mapping.originInfo());
  assert(source.vectorSize() == mapping.restrictedInfo());
  
  Vector buffer(mapping.originInfo());
  for (VecBasis::const_iterator vecIt = source.begin(); vecIt != source.end(); ++vecIt) {
    mapping.extension(*vecIt, buffer);
    out << buffer;
  }
  return out;
}

BasisOutputStream &
writeRestrictedVectors(BasisOutputStream &out, const VecBasis &source, const NodalRestrictionMapping &mapping) {
  assert(out.vectorSize() == mapping.restrictedInfo());
  assert(source.vectorSize() == mapping.originInfo());
  
  Vector buffer(mapping.restrictedInfo());
  for (VecBasis::const_iterator vecIt = source.begin(); vecIt != source.end(); ++vecIt) {
    mapping.restriction(*vecIt, buffer);
    out << buffer;
  }
  return out;
}

// Input
BasisInputStream &
operator>>(BasisInputStream &in, VecBasis &sink) {
  sink.dimensionIs(in.size() - in.currentVectorRank(), in.vectorSize());
  return readVectors(in, sink.begin(), sink.end());
}

BasisInputStream &
readVectors(BasisInputStream &in, VecBasis &sink, int countMax) {
  const int count = std::max(std::min(in.size() - in.currentVectorRank(), countMax), 0);
  sink.dimensionIs(count, in.vectorSize());
  return readVectors(in, sink.begin(), sink.end());
}

BasisInputStream &
readRestrictedVectors(BasisInputStream &in, VecBasis &sink, const NodalRestrictionMapping &mapping) {
  assert(in.vectorSize() == mapping.originInfo());
  
  sink.dimensionIs(in.size() - in.currentVectorRank(), mapping.restrictedInfo());
  Vector buffer(mapping.originInfo());
  for (VecBasis::iterator vecIt = sink.begin(); vecIt != sink.end(); ++vecIt) {
    in >> buffer;
    mapping.restriction(buffer, *vecIt);
  }
  return in;
}

BasisInputStream &
readExtendedVectors(BasisInputStream &in, VecBasis &sink, const NodalRestrictionMapping &mapping) {
  assert(in.vectorSize() == mapping.restrictedInfo());
  
  sink.dimensionIs(in.size() - in.currentVectorRank(), mapping.originInfo());
  Vector buffer(mapping.restrictedInfo());
  for (VecBasis::iterator vecIt = sink.begin(); vecIt != sink.end(); ++vecIt) {
    in >> buffer;
    mapping.extension(buffer, *vecIt);
  }
  return in;
}

} /* end namespace Rom */
