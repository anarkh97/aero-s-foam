#ifndef ROM_VECBASIS_H
#define ROM_VECBASIS_H

#include <memory>
#include <algorithm>
#include <Feti.d/DistrVector.h>
#include <Threads.d/PHelper.h>
#include <Utils.d/dofset.h>
#ifdef USE_EIGEN3
#include <Eigen/Core>
#include <Eigen/Sparse>
#endif

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
#ifdef USE_EIGEN3
  Eigen::Matrix<Scalar,Eigen::Dynamic,Eigen::Dynamic, Eigen::RowMajor> V;
#endif
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

  //GenVecType<Scalar> & project(GenVecType<Scalar> &, GenVecType<Scalar> &);
  //void domainMatrixVecMult(GenVecType<Scalar> &, GenVecType<Scalar> &);
/*
  GenDistrVector<double> & project(GenDistrVector<double> &, GenDistrVector<double> &);
  GenDistrVector<double> & projectUp(GenDistrVector<double> &, GenDistrVector<double> &);
  GenDistrVector<double> & projectUp(std::vector<double> &   , GenDistrVector<double> & );
  GenDistrVector<double> & projectDown(GenDistrVector<double> &, GenDistrVector<double> &);
*/
  VecType & project(VecType &, VecType &) const;
  VecType & projectUp(VecType &, VecType &) const;
  VecType & projectUp(std::vector<double> &, VecType &) const;
  VecType & projectDown(VecType &, VecType &) const;
 
  void makeSparseBasis(std::vector<int> &, DofSetArray *); 

  timespec tS1, tS2;
  double time1, time2, time3, time4, time5, time6, counter;

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
#ifdef USE_EIGEN3
  Eigen::Map<Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> > basis;
  Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> compressedBasis; 
  std::vector<int> compressedKey;
#endif
};

template <typename Scalar, template <typename Scalar> class GenVecType>
inline
void
GenVecBasis<Scalar, GenVecType>::swap(GenVecBasis &other) {
  std::swap(vectorInfo_,  other.vectorInfo_);
  std::swap(vectorCount_, other.vectorCount_);
  std::swap(buffer_,      other.buffer_);
  std::swap(vectors_,     other.vectors_);

#ifdef USE_EIGEN3
  new (&basis) Eigen::Map< Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> >(buffer_, vectorSize(), vectorCount());
#endif

  time1 = 0; time2 = 0; time3 = 0; time4 = 0; time5 = 0; time6 = 0; counter = 0;
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
#ifdef USE_EIGEN3
  new (&basis) Eigen::Map< Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> >(buffer_, vectorSize(), vectorCount());
#endif
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
#ifdef USE_EIGEN3
 ,basis(NULL,0,0),
 compressedBasis(0,0),
 compressedKey(0)
#endif
{
  placeVectors();
}



template <typename Scalar, template <typename Scalar> class GenVecType>
GenVecBasis<Scalar, GenVecType>::GenVecBasis(int vCount, InfoType vInfo) :
 vectorInfo_(vInfo),
 vectorCount_(vCount)
#ifdef USE_EIGEN3
 ,basis(NULL,0,0),
 compressedBasis(0,0),
 compressedKey(0)
#endif
{
  placeVectors();
}

template <typename Scalar, template <typename Scalar> class GenVecType>
GenVecBasis<Scalar, GenVecType>::GenVecBasis(const GenVecBasis &other) :
 vectorInfo_(other.vectorInfo_),
 vectorCount_(other.vectorCount_)
#ifdef USE_EIGEN3
 ,basis(NULL,0,0),
 compressedBasis(0,0),
 compressedKey(0)
#endif
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
