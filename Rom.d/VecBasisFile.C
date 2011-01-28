#include "VecBasisFile.h"

#include "BasisFileStream.h"

#include <Math.d/Vector.h>

BasisInputStream &
operator>>(BasisInputStream &in, VecBasis &sink) {
  GenVecBasis<double> temp(in.size(), in.vectorSize());

  GenVecBasis<double>::iterator itEnd = temp.end();
  for (GenVecBasis<double>::iterator it = temp.begin(); it != itEnd; ++it) {
    in >> *it;
  }
  sink.swap(temp);

  return in;
}

BasisOutputStream &
operator<<(BasisOutputStream &out, const VecBasis &source) {
  VecBasis::const_iterator itEnd = source.end();
  for (VecBasis::const_iterator it = source.begin(); it != itEnd; ++it) {
    out << *it;
  }
  return out;
}
