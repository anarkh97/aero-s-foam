#ifndef ROM_BASISFILEITERATOR_H
#define ROM_BASISFILEITERATOR_H

#include "BasisInputFile.h"
#include "BasisOutputFile.h"
#include "VecNodeDof6Conversion.h"

#include "NodeDof6Buffer.h"

#include <string>
#include <utility>
#include <cstddef>

namespace Rom {

// Input Stream

class BasisInputStream {
public:
  int size() const { return file_.stateCount(); }
  int vectorSize() const { return converter_.vectorSize(); }
  int nodeCount() const { return converter_.nodeCount(); }

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
  template <typename VectorBufferType>
  void performInput(VectorBufferType &);
  
  template <typename VectorBufferType>
  void performInput(std::pair<double, VectorBufferType> &);

  BasisInputFile file_;
  const VecNodeDof6Conversion &converter_;
  NodeDof6Buffer buffer_;
};

inline
BasisInputStream::operator const void*() const {
  return file_.validCurrentState() ? this : NULL;
}

template <typename VectorBufferType>
inline
void
BasisInputStream::performInput(VectorBufferType &target) {
  converter_.vector(file_.currentStateBuffer(buffer_), target);
  file_.currentStateIndexInc();
}
 
template <typename VectorBufferType>
inline
void
BasisInputStream::performInput(std::pair<double, VectorBufferType> &target) {
  target.first = file_.currentStateHeaderValue();
  performInput(target.second);
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
  while (it != last && in) {
    in >> *it++;
  }
  return in;
}

// Output Stream

class BasisOutputStream {
public:
  int size() const { return file_.stateCount(); }
  int vectorSize() const { return converter_.vectorSize(); }
  int nodeCount() const { return converter_.nodeCount(); }
  
  BasisOutputStream(const std::string &fileName, const VecNodeDof6Conversion &converter);

  template <typename VectorBufferType>
  friend BasisOutputStream &operator<<(BasisOutputStream &, const VectorBufferType &source);
  
  template <typename VectorBufferType>
  friend BasisOutputStream &operator<<(BasisOutputStream &, const std::pair<double, VectorBufferType> &source);
  
private:
  template <typename VectorBufferType>
  const NodeDof6Buffer &convert(const VectorBufferType &);

  BasisOutputFile file_;
  const VecNodeDof6Conversion &converter_;

  NodeDof6Buffer buffer_;
};

// Output operations

template <typename VectorBufferType>
inline
const NodeDof6Buffer &
BasisOutputStream::convert(const VectorBufferType &source) {
  return converter_.nodeDof6(source, buffer_);
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
