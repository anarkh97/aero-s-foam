#ifndef ROM_BASISFILEITERATOR_H
#define ROM_BASISFILEITERATOR_H

#include "BasisInputFile.h"
#include "BasisOutputFile.h"
#include "VecNodeDof6Conversion.h"

#include "NodeDof6Buffer.h"

#include <Math.d/Vector.h>

#include <string>
#include <utility>
#include <cstddef>

// Input Iterator

class BasisInputRange {
public:
  int size() const { return file_.stateCount(); }
  int vectorSize() const { return converter_.vectorSize(); }
  int nodeCount() const { return converter_.nodeCount(); }

  class const_iterator;
  const_iterator begin();
  const_iterator end() const;

  BasisInputRange(const std::string &fileName, const VecNodeDof6Conversion &converter);

  class Copyer;

private:
  BasisInputFile file_;
  const VecNodeDof6Conversion &converter_;

  NodeDof6Buffer buffer_;

  friend class const_iterator;
  friend class Copyer;
};

class BasisInputRange::Copyer {
public:
  template <typename BufferType>
  void operator()(const BufferType &targetBuffer) const {
    parent_.converter_.vector(parent_.buffer_, targetBuffer);
  }
  
private:
  explicit Copyer(const BasisInputRange &parent) :
    parent_(parent)
  {}

  const BasisInputRange &parent_;

  friend class BasisInputRange::const_iterator;
};

class BasisInputRange::const_iterator {
public:
  BasisInputRange::Copyer operator*() const {
    return BasisInputRange::Copyer(*parent_);
  }

  const_iterator &operator++();

  bool operator==(const const_iterator &other) {
    return parent_ == other.parent_;
  }

  bool operator!=(const const_iterator &other) {
    return !(*this == other);
  }

  const_iterator() :
    parent_(NULL)
  {} 

private:
  explicit const_iterator(BasisInputRange *parent) :
    parent_(parent)
  {
    validateParent();
  }

  void validateParent();

  BasisInputRange *parent_;

  friend class BasisInputRange;
};

inline
BasisInputRange::const_iterator
BasisInputRange::begin() {
  return const_iterator(this);
}

inline
BasisInputRange::const_iterator
BasisInputRange::end() const {
  return const_iterator();
}

// Output Iterator

class BasisOutputRange {
public:
  int size() const { return file_.stateCount(); }
  int vectorSize() const { return converter_.vectorSize(); }
  int nodeCount() const { return converter_.nodeCount(); }
  
  class iterator;
  iterator begin();
  iterator end() const;

  BasisOutputRange(const std::string &fileName, const VecNodeDof6Conversion &converter);

  class Copyer;

private:
  BasisOutputFile file_;
  const VecNodeDof6Conversion &converter_;

  NodeDof6Buffer buffer_;

  friend class iterator;
  friend class Copyer;
}; 

class BasisOutputRange::Copyer {
public:
  template <typename BufferType>
  void operator=(const std::pair<double, BufferType> &entry);

private:
  explicit Copyer(BasisOutputRange &parent) :
    parent_(parent)
  {}

  BasisOutputRange &parent_;

  friend class BasisOutputRange::iterator;
};

template <typename BufferType>
void
BasisOutputRange::Copyer::operator=(const std::pair<double, BufferType> &entry) {
  const NodeDof6Buffer &buffer = parent_.converter_.nodeDof6(entry.second, parent_.buffer_);
  parent_.file_.stateAdd(buffer, entry.first);
}

class BasisOutputRange::iterator {
public:
  BasisOutputRange::Copyer operator*() {
    return BasisOutputRange::Copyer(*parent_);
  }

  iterator &operator++() {
    // Nothing to do 
    return *this;
  }

  iterator() :
    parent_(NULL)
  {}

private:
  explicit iterator(BasisOutputRange &parent) :
    parent_(&parent)
  {}

  BasisOutputRange *parent_;

  friend class BasisOutputRange;
};

inline
BasisOutputRange::iterator
BasisOutputRange::begin() {
  return iterator(*this);
}

inline
BasisOutputRange::iterator
BasisOutputRange::end() const {
  return iterator();
}

#endif /* ROM_BASISFILEITERATOR_H */
