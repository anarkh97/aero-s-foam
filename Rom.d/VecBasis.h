#ifndef ROM_VECBASIS_H
#define ROM_VECBASIS_H

template <typename Scalar> class GenVector;

#include <memory>

template <typename Scalar, template <typename Scalar> class GenVecType = GenVector>
class GenVecBasis : private std::allocator<GenVecType<Scalar> > {
public:
  typedef GenVecType<Scalar> VecType;

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

  // Ctor and dtor
  GenVecBasis(int vCount, int vSize);
  ~GenVecBasis();

private:
  typedef std::allocator<VecType> Allocator;

  int vectorSize_;
  int vectorCount_;

  Scalar *buffer_;
  VecType *vectors_;

  // Disallow copy and assignment
  GenVecBasis(const GenVecBasis &);
  GenVecBasis &operator=(const GenVecBasis &);
};

template <typename Scalar, template <typename Scalar> class GenVecType>
GenVecBasis<Scalar, GenVecType>::GenVecBasis(int vCount, int vSize) :
 vectorSize_(vSize),
 vectorCount_(vCount),
 buffer_(new Scalar[vectorSize_ * vectorCount_]),
 vectors_(Allocator::allocate(vectorCount_))
{
  for (int iVec = 0; iVec < vectorCount_; ++iVec) {
    new(vectors_ + iVec) VecType(vectorSize_, buffer_ + (iVec * vectorSize_), false);
  }
}

template <typename Scalar, template <typename Scalar> class GenVecType>
GenVecBasis<Scalar, GenVecType>::~GenVecBasis() {
  for (int iVec = 0; iVec < vectorCount_; ++iVec) {
    Allocator::destroy(vectors_ + iVec);
  }

  Allocator::deallocate(vectors_, vectorCount_);
  delete[] buffer_;
}
  
typedef GenVecBasis<double> VecBasis;

#endif /* ROM_VECBASIS_H */
