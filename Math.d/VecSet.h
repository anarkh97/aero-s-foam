#ifndef _VECTOR_SET_H_
#define _VECTOR_SET_H_

#ifdef OLD_STL
#include <memory.h>
#ifdef USE_IOSTREAM
#include <iostream>
using std::cerr;
using std::endl;
#endif
#else
#include <memory>
using std::allocator;
#ifdef USE_IOSTREAM
#include <iostream>
using std::cerr;
using std::endl;
#endif
#endif

//------------------------------------------------------------------------------

template<class VecType>
class VecSet {

  int numVec;
  const typename VecType::InfoType &len;

  allocator<VecType> alloc;
  VecType *vecSet;

public:

  VecSet(int, const typename VecType::InfoType &);
  VecSet(const VecSet<VecType> &);
  ~VecSet();

  VecType &operator[] (int i) const { return vecSet[i]; }
  
  int numVectors() const { return numVec; }
  typename VecType::InfoType &size() const { return len; }

  void resize(int);

#ifdef USE_IOSTREAM
  void print(char *msg = "") {
    if (msg) cerr << msg << endl;

    for (int i=0; i<numVec; ++i) {
      cerr << "vector " << i << ":";
      vecSet[i].print();
    }
  }
#endif

};

//------------------------------------------------------------------------------

template<class VecType>
VecSet<VecType>::VecSet(int _numVec, const typename VecType::InfoType &_len) : len(_len)
{

  numVec = _numVec;

  vecSet = static_cast<VecType*>(alloc.allocate(numVec));

  for (int i = 0; i < numVec; ++i)
    new (static_cast<void *>(vecSet+i)) VecType(len);

}

//------------------------------------------------------------------------------

template<class VecType>
void VecSet<VecType>::resize(int n)
{
 int i;
 for (i = 0; i < numVec; ++i)
   vecSet[i].~VecType();

 alloc.deallocate(vecSet, numVec);

 numVec = n;

 vecSet = static_cast<VecType*>(alloc.allocate(numVec));

 for (i = 0; i < numVec; ++i)
    new (static_cast<void *>(vecSet+i)) VecType(len);

}

//------------------------------------------------------------------------------

template<class VecType>
VecSet<VecType>::VecSet(const VecSet<VecType> &vectorSet)
{
 
  numVec = vectorSet.numVectors();
  len = vectorSet.size();

  vecSet = static_cast<VecType*>(alloc.allocate(numVec));

  for (int i = 0; i < numVec; ++i) {
    new (static_cast<void *>(vecSet+i)) VecType(len);
    vecSet[i] = vectorSet[i];
  }

}

//------------------------------------------------------------------------------

template<class VecType>
VecSet<VecType>::~VecSet() 
{
  if (vecSet) {
    for (int i = 0; i < numVec; ++i) vecSet[i].~VecType();
#ifdef OLD_STL
    alloc.deallocate(vecSet);
#else
    alloc.deallocate(vecSet, numVec);
#endif
    vecSet = 0;
  }
}

//------------------------------------------------------------------------------

#endif
