#ifndef ROM_VECBASIS_H
#define ROM_VECBASIS_H

#include <memory>
#include <algorithm>
#include <Feti.d/DistrVector.h>
#include <Threads.d/PHelper.h>
#include <Eigen/Core>
#include <Eigen/Sparse>

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
  Eigen::Matrix<Scalar,Eigen::Dynamic,Eigen::Dynamic, Eigen::RowMajor> V;
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

  GenVecType<Scalar> & project(GenVecType<Scalar> &, GenVecType<Scalar> &);
  void domainMatrixVecMult(GenVecType<Scalar> &, GenVecType<Scalar> &);

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
GenVecType<Scalar> &
GenVecBasis<Scalar, GenVecType>::project(GenVecType<Scalar> &x, GenVecType<Scalar> &result) {

   std::cout << " ... model3 ..." << std::endl;
   domainMatrixVecMult( x, result);

   return result;
}

template <typename Scalar, template <typename Scalar> class GenVecType>
inline
void
GenVecBasis<Scalar, GenVecType>::domainMatrixVecMult(GenVecType<Scalar> &_f, GenVecType<Scalar> &_result) {

  Eigen::Matrix<Scalar, Eigen::Dynamic, 1> GenCoordinates; 
  Eigen::Map< Eigen::Matrix<Scalar, Eigen::Dynamic, 1> > result(_result.data(), _result.size());
  Eigen::Map< Eigen::Matrix<Scalar, Eigen::Dynamic, 1> > f(_f.data(), size());

  typedef Eigen::Triplet<Scalar> T;
  std::vector<T> tripletList;

  std::cout << "line 171" << std::endl;
  for (int i = 0; i < size(); i++) {
   if(f(i) != 0.){
     std::cout << "element " << i << " = " << f(i) << std::endl;
     tripletList.push_back(T(i,1,f(i)));}
  }
  std::cout << "list size = " << tripletList.size() << std::endl;
  std::cout << " line 178 " << std::endl;
  Eigen::SparseMatrix<Scalar> sparsef(size(),1);
  std::cout << " line 180 " << std::endl;
  sparsef.setFromTriplets(tripletList.begin(), tripletList.end());
  std::cout << "line 182 " << std::endl;
  Eigen::Map< Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic>, Eigen::Unaligned, Eigen::OuterStride<> > basisT(vectors_[0].data(), vectors_[0].size(), vectorCount(), Eigen::OuterStride<>(size()));

  std::cout << "line 184" << std::endl;
  GenCoordinates = basisT*sparsef;
  std::cout << "line 186" << std::endl;

  if(structCom)
    structCom->globalSum(GenCoordinates.size(), GenCoordinates.data());

  result = basisT.transpose()*GenCoordinates;

}

/*template <typename Scalar, template <typename Scalar> class GenVecType>
void
GenVecBasis<Scalar, GenVecType>::domainMatrixVecMultPar(GenVecType<Scalar> &f, GenVector<Scalar> &result) {

  GenVector<Scalar> *subresult = new GenVector<Scalar>[f.num()];
  execParal2R(f.num() , const_cast<GenVecBasis<Scalar, GenVecType> *>(this), &GenVecBasis<Scalar, GenVecType>::subdomainMatrixVecMult, f, subresult);

  result.zero();
  for(int i=0; i<numSub; ++i)
    result += subresult[i];
  if(structCom)
      structCom->globalSum(result.size(), result.data());

  delete [] subresult;
}

template <typename Scalar, template <typename Scalar> class GenVecType>
void
GenVecBasis<Scalar, GenVecType>::subdomainMatrixVecMult(int isub, GenVecType<Scalar> &f, GenVector<Scalar> *result) {

  Eigen::Map< Eigen::Matrix<Scalar, Eigen::Dynamic, 1> > RES(result[isub].data(), result[isub].size()); 
  Eigen::Map< Eigen::Matrix<Scalar, Eigen::Dynamic, 1> > fSub(f.subData(iSub), size());

  Eigen::Map< Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic>, Eigen::Unaligned, Eigen::OuterStride<> > basisSubBlock(vectors_[0].subData(iSub), vectors_[0].subLen(iSub), vectorCount(), Eigen::OuterStride<>(size()));

  RES = basisSubBlock*fSub;
}*/

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
