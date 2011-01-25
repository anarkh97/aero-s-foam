#include "BasisFileIterator.h"

#include <stdexcept>

// Input Iterator

BasisInputRange::BasisInputRange(const std::string &fileName, const VecNodeDof6Conversion &converter) :
  file_(fileName),
  converter_(converter),
  buffer_(file_.nodeCount())
{
  if (file_.nodeCount() != converter_.nodeCount()) {
    throw std::invalid_argument("Incompatible node counts");
  }
}

BasisInputRange::const_iterator &
BasisInputRange::const_iterator::operator++() {
  parent_->file_.currentStateIndexInc();
  validateParent();

  return *this;
}

void
BasisInputRange::const_iterator::validateParent() {
  if (parent_->file_.validCurrentState()) {
    parent_->file_.currentStateBuffer(parent_->buffer_);
  } else {
    parent_ = NULL;
  }
}

// Output Iterator

BasisOutputRange::BasisOutputRange(const std::string &fileName, const VecNodeDof6Conversion &converter) :
  file_(fileName, converter.nodeCount()),
  converter_(converter),
  buffer_(file_.nodeCount())
{}

