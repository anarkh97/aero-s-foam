#include "BasisFileStream.h"

#include <stdexcept>

namespace Rom {

BasisInputStream::BasisInputStream(const std::string &fileName, const VecNodeDof6Conversion &converter) :
  file_(fileName),
  isValid_(true),
  converter_(converter),
  buffer_(file_.nodeIdBegin(), file_.nodeIdEnd())
{}


BasisOutputStream::BasisOutputStream(const std::string &fileName, const VecNodeDof6Conversion &converter, bool restart) :
  file_(fileName, converter.dofSetNodeCount(), restart), // Conservative, potentially overallocating
  converter_(converter),
  buffer_(file_.nodeCount())
{}

} /* end namespace Rom */
