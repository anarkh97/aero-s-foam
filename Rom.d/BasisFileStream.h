#ifndef ROM_BASISFILEITERATOR_H
#define ROM_BASISFILEITERATOR_H

#include "BasisBinaryFile.h"
#include "VecNodeDof6Conversion.h"

#include "NodeDof6Buffer.h"
#include "MappedNodeDof6Buffer.h"

#include <string>
#include <utility>
#include <cstddef>

namespace Rom {

// Input Stream

class BasisInputStream {
public:
  int size() const { return file_.stateCount(); }
  int vectorSize() const { return converter_.vectorSize(); }
  int dofSetNodeCount() const { return converter_.dofSetNodeCount(); }

  // Check the readiness of the stream
  int currentVectorRank() const { return file_.currentStateIndex(); }
  operator const void*() const;

  BasisInputStream(const std::string &fileName, const VecNodeDof6Conversion &converter);

  template <typename VectorBufferType>
  friend BasisInputStream &operator>>(BasisInputStream &, VectorBufferType &);
  
  // Convenience overload (to bind with non-const rvalues)
  friend BasisInputStream &operator>>(BasisInputStream &, double *);
  
  template <typename VectorBufferType>
  friend BasisInputStream &operator>>(BasisInputStream &, std::pair<double, VectorBufferType> &);

private:
  template <typename VectorBufferType> void performInput(VectorBufferType &);
  template <typename VectorBufferType> void performUncheckedInput(VectorBufferType &);
 
  void checkInput();

  template <typename VectorBufferType>
  void performInput(std::pair<double, VectorBufferType> &);

  BasisBinaryInputFile file_;
  bool isValid_;
  const VecNodeDof6Conversion &converter_;
  MappedNodeDof6Buffer buffer_;
};

inline
BasisInputStream::operator const void*() const {
  return isValid_ ? this : NULL;
}

inline
void
BasisInputStream::checkInput() {
  isValid_ = file_.validCurrentState();
}

template <typename VectorBufferType>
inline
void
BasisInputStream::performUncheckedInput(VectorBufferType &target) {
  file_.currentStateBuffer(buffer_.underlyingBuffer()); // Fill buffer in one pass, ignoring the file node mapping
  converter_.vector(buffer_, target); // Fill vector according to the file node mapping 
  file_.currentStateIndexInc();
}

template <typename VectorBufferType>
inline
void
BasisInputStream::performInput(VectorBufferType &target) {
  checkInput();
  if (isValid_) {
    performUncheckedInput(target);
  }
}
 
template <typename VectorBufferType>
inline
void
BasisInputStream::performInput(std::pair<double, VectorBufferType> &target) {
  checkInput();
  if (isValid_) {
    target.first = file_.currentStateHeaderValue();
    performUncheckedInput(target.second);
  }
}

// Input operations

template <typename VectorBufferType>
BasisInputStream &
operator>>(BasisInputStream &in, VectorBufferType &target) {
  in.performInput(target);
  return in;
}

inline
BasisInputStream &
operator>>(BasisInputStream &in, double *target) {
  in.performInput(target);
  return in;
}

template <typename VectorBufferType>
BasisInputStream &
operator>>(BasisInputStream &in, std::pair<double, VectorBufferType> &target) {
  in.performInput(target);
  return in;
}

template <typename FwdIter>
BasisInputStream &
readVectors(BasisInputStream &in, FwdIter first, FwdIter last) {
  FwdIter it = first;
  while (it != last && in >> *it++) {
    // Nothing to do
  }
  return in;
}

// Output Stream

class BasisOutputStream {
public:
  int size() const { return file_.stateCount(); }
  int vectorSize() const { return converter_.vectorSize(); }
  int dofSetNodeCount() const { return converter_.dofSetNodeCount(); }
  
  BasisOutputStream(const std::string &fileName, const VecNodeDof6Conversion &converter, bool restart);

  template <typename VectorBufferType>
  friend BasisOutputStream &operator<<(BasisOutputStream &, const VectorBufferType &source);
  
  template <typename VectorBufferType>
  friend BasisOutputStream &operator<<(BasisOutputStream &, const std::pair<double, VectorBufferType> &source);
  
private:
  template <typename VectorBufferType>
  const NodeDof6Buffer &convert(const VectorBufferType &);

  BasisBinaryOutputFile file_;
  const VecNodeDof6Conversion &converter_;

  NodeDof6Buffer buffer_;
};

// Output operations

template <typename VectorBufferType>
inline
const NodeDof6Buffer &
BasisOutputStream::convert(const VectorBufferType &source) {
  return converter_.paddedNodeDof6(source, buffer_);
}

template <typename VectorBufferType>
BasisOutputStream &
operator<<(BasisOutputStream &out, const std::pair<double, VectorBufferType> &source) {
  out.file_.stateAdd(out.convert(source.second), source.first);
  return out;
}

template <typename VectorBufferType>
BasisOutputStream &
operator<<(BasisOutputStream &out, const VectorBufferType &source) {
  out.file_.stateAdd(out.convert(source));
  return out;
}

template <typename InputIter>
BasisOutputStream &
writeVectors(BasisOutputStream &out, InputIter first, InputIter last) {
  InputIter it = first;
  while (it != last) {
    out << *it++;
  }
  return out;
}

} /* end namespace Rom */

#endif /* ROM_BASISFILEITERATOR_H */
