#ifndef ROM_VECBASIS_H
#define ROM_VECBASIS_H

#include <memory>
#include <algorithm>

template <typename Scalar> class GenVector;

namespace Rom {

template <typename Scalar, template <typename Scalar> class GenVecType = GenVector>
class GenVecBasis : private std::allocator<GenVecType<Scalar> > {
public:
  typedef GenVecType<Scalar> VecType;
  typedef typename VecType::InfoType InfoType;

  // Immutable individual vector size (and compatibility aliases)
  int vectorSize() const { return vectorSize_; }
  int size()       const { return vectorSize_; }

  // Immutable vector count (and compatibility aliases)
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
  GenVecBasis(int vCount, int vSize);
 
  // Copy, assignment and swap 
  GenVecBasis(const GenVecBasis &);
  GenVecBasis &operator=(const GenVecBasis &);
  void swap(GenVecBasis &);
  
  // Reshaping
  void dimensionIs(int vCount, int vSize);

  ~GenVecBasis();

private:
  typedef std::allocator<VecType> Allocator;

  void placeVectors();
  void copyBufferContent(const GenVecBasis &other);

  int vectorSize_;
  int vectorCount_;

  Scalar *buffer_;
  VecType *vectors_;
};

template <typename Scalar, template <typename Scalar> class GenVecType>
inline
void
GenVecBasis<Scalar, GenVecType>::swap(GenVecBasis &other) {
  std::swap(vectorSize_,  other.vectorSize_);
  std::swap(vectorCount_, other.vectorCount_);
  std::swap(buffer_,      other.buffer_);
  std::swap(vectors_,     other.vectors_);
}

template <typename Scalar, template <typename Scalar> class GenVecType>
inline
void
GenVecBasis<Scalar, GenVecType>::placeVectors() {
  buffer_ = new Scalar[vectorSize_ * vectorCount_];
  vectors_ = Allocator::allocate(vectorCount_);
  for (int iVec = 0; iVec < vectorCount_; ++iVec) {
    new(vectors_ + iVec) VecType(vectorSize_, buffer_ + (iVec * vectorSize_), false);
  }
}

template <typename Scalar, template <typename Scalar> class GenVecType>
inline
void
GenVecBasis<Scalar, GenVecType>::copyBufferContent(const GenVecBasis &other) {
  std::copy(other.buffer_, other.buffer_ + vectorSize_ * vectorCount_, buffer_);
}

template <typename Scalar, template <typename Scalar> class GenVecType>
GenVecBasis<Scalar, GenVecType>::GenVecBasis() :
 vectorSize_(),
 vectorCount_(0)
{
  placeVectors();
}

template <typename Scalar, template <typename Scalar> class GenVecType>
GenVecBasis<Scalar, GenVecType>::GenVecBasis(int vCount, int vSize) :
 vectorSize_(vSize),
 vectorCount_(vCount)
{
  placeVectors();
}

template <typename Scalar, template <typename Scalar> class GenVecType>
GenVecBasis<Scalar, GenVecType>::GenVecBasis(const GenVecBasis &other) :
 vectorSize_(other.vectorSize_),
 vectorCount_(other.vectorCount_)
{
  placeVectors();
  copyBufferContent(other);
}

template <typename Scalar, template <typename Scalar> class GenVecType>
GenVecBasis<Scalar, GenVecType> &
GenVecBasis<Scalar, GenVecType>::operator=(const GenVecBasis &other) {
  if (this != &other) {
    if (vectorCount_ == other.vectorCount_ && vectorSize_ == other.vectorSize_) {
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
GenVecBasis<Scalar, GenVecType>::dimensionIs(int vCount, int vSize) {
  if (vCount != vectorCount() || vSize != vectorSize()) {
    GenVecBasis temp(vCount, vSize);
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
