#include "BasisFileStream.h"

#include <stdexcept>

namespace Rom {

// Input Iterator

BasisInputStream::BasisInputStream(const std::string &fileName, const VecNodeDof6Conversion &converter) :
  file_(fileName),
  isValid_(true),
  converter_(converter),
  buffer_(file_.nodeCount())
{
  if (file_.nodeCount() != converter_.nodeCount()) {
    throw std::invalid_argument("Incompatible node counts");
  }
}

// Output Iterator

BasisOutputStream::BasisOutputStream(const std::string &fileName, const VecNodeDof6Conversion &converter) :
  file_(fileName, converter.nodeCount()),
  converter_(converter),
  buffer_(file_.nodeCount())
{}

} /* end namespace Rom */
