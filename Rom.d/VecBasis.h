#ifndef ROM_VECBASIS_H
#define ROM_VECBASIS_H

#include <memory>
#include <algorithm>

template <typename Scalar> class GenVector;

namespace Rom {

template <typename Scalar, template <typename Scalar> class GenVecType>
struct VecTraits {
  typedef GenVecType<Scalar> Type;
  typedef typename Type::InfoType InfoType;
  typedef InfoType InternalInfoType;
  
  static InfoType defaultInfo() { return InfoType(); }
  static int length(InfoType info) { return info; }
  static bool equals(InfoType i, InfoType j) { return i == j; }
  static bool not_equals(InfoType i, InfoType j) { return i != j; }
};

template <typename Scalar, template <typename Scalar> class GenVecType = GenVector>
class GenVecBasis : private std::allocator<GenVecType<Scalar> > {
private:
  typedef VecTraits<Scalar, GenVecType> Traits;
public:
  typedef typename Traits::Type VecType;
  typedef typename Traits::InfoType InfoType;

  // Common vector information
  InfoType vectorInfo() const { return vectorInfo_; }

  // Individual vector size (and compatibility aliases)
  int vectorSize() const { return Traits::length(vectorInfo_); }
  int size()       const { return vectorSize(); }

  // Vector count (and compatibility aliases)
  int vectorCount() const { return vectorCount_; }
  int numVec()      const { return vectorCount_; }
  int numVectors()  const { return vectorCount_; }

  // Iteration
  typedef const VecType *const_iterator;
  const_iterator begin() const { return vectors_; }
  const_iterator end() const { return vectors_ + vectorCount_; }
  
  typedef VecType *iterator;
  iterator begin() { return vectors_; }
  iterator end() { return vectors_ + vectorCount_; }

  // Unchecked direct individual vector read access
  const VecType &operator[](int i) const { return vectors_[i]; }

  // Unchecked direct individual vector write access
  // (must take care to NOT reallocate underlying memory)
  VecType &operator[](int i) { return vectors_[i]; }

  // Constructors
  GenVecBasis();
  GenVecBasis(int vCount, InfoType vInfo);
 
  // Copy, assignment and swap 
  GenVecBasis(const GenVecBasis &);
  GenVecBasis &operator=(const GenVecBasis &);
  void swap(GenVecBasis &);
  
  // Reshaping
  void dimensionIs(int vCount, InfoType vInfo);

  ~GenVecBasis();

private:
  typedef std::allocator<VecType> Allocator;

  void placeVectors();
  void copyBufferContent(const GenVecBasis &other);

  typename Traits::InternalInfoType vectorInfo_;
  int vectorCount_;

  Scalar *buffer_;
  VecType *vectors_;
};

template <typename Scalar, template <typename Scalar> class GenVecType>
inline
void
GenVecBasis<Scalar, GenVecType>::swap(GenVecBasis &other) {
  std::swap(vectorInfo_,  other.vectorInfo_);
  std::swap(vectorCount_, other.vectorCount_);
  std::swap(buffer_,      other.buffer_);
  std::swap(vectors_,     other.vectors_);
}

template <typename Scalar, template <typename Scalar> class GenVecType>
inline
void
GenVecBasis<Scalar, GenVecType>::placeVectors() {
  buffer_ = new Scalar[vectorSize() * vectorCount_];
  vectors_ = Allocator::allocate(vectorCount_);
  for (int iVec = 0; iVec < vectorCount_; ++iVec) {
    new(vectors_ + iVec) VecType(vectorInfo_, buffer_ + (iVec * vectorSize()), false);
  }
}

template <typename Scalar, template <typename Scalar> class GenVecType>
inline
void
GenVecBasis<Scalar, GenVecType>::copyBufferContent(const GenVecBasis &other) {
  std::copy(other.buffer_, other.buffer_ + vectorSize() * vectorCount_, buffer_);
}

template <typename Scalar, template <typename Scalar> class GenVecType>
GenVecBasis<Scalar, GenVecType>::GenVecBasis() :
 vectorInfo_(Traits::defaultInfo()),
 vectorCount_(0)
{
  placeVectors();
}

template <typename Scalar, template <typename Scalar> class GenVecType>
GenVecBasis<Scalar, GenVecType>::GenVecBasis(int vCount, InfoType vInfo) :
 vectorInfo_(vInfo),
 vectorCount_(vCount)
{
  placeVectors();
}

template <typename Scalar, template <typename Scalar> class GenVecType>
GenVecBasis<Scalar, GenVecType>::GenVecBasis(const GenVecBasis &other) :
 vectorInfo_(other.vectorInfo_),
 vectorCount_(other.vectorCount_)
{
  placeVectors();
  copyBufferContent(other);
}

template <typename Scalar, template <typename Scalar> class GenVecType>
GenVecBasis<Scalar, GenVecType> &
GenVecBasis<Scalar, GenVecType>::operator=(const GenVecBasis &other) {
  if (this != &other) {
    if (vectorCount_ == other.vectorCount_ && Traits::equals(vectorInfo_, other.vectorInfo_)) {
      copyBufferContent(other);
    } else {
      GenVecBasis temp(other);
      swap(temp);
    }
  }
  return *this;
}

template <typename Scalar, template <typename Scalar> class GenVecType>
void
GenVecBasis<Scalar, GenVecType>::dimensionIs(int vCount, InfoType vInfo) {
  if (vCount != vectorCount_ || Traits::not_equals(vInfo, vectorInfo_)) {
    GenVecBasis temp(vCount, vInfo);
    swap(temp);
  }
}

template <typename Scalar, template <typename Scalar> class GenVecType>
GenVecBasis<Scalar, GenVecType>::~GenVecBasis() {
  int iVec = vectorCount_;
  while (iVec--) {
    Allocator::destroy(vectors_ + iVec);
  }

  Allocator::deallocate(vectors_, vectorCount_);
  delete[] buffer_;
}

template <typename Scalar, template <typename Scalar> class GenVecType>
inline
void
swap(GenVecBasis<Scalar, GenVecType> &a, GenVecBasis<Scalar, GenVecType> &b) {
  a.swap(b);
}

typedef GenVecBasis<double> VecBasis;

} /* end namespace Rom */

#endif /* ROM_VECBASIS_H */
