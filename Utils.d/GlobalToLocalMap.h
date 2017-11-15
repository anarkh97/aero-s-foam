#ifndef _GLOBALTOLOCALMAP_H_
#define _GLOBALTOLOCALMAP_H_

//#define HB_USE_MEMCPY
#ifdef HB_USE_MEMCPY
#include <cstring>
#endif

//#define MAP_MIN_MEMORY
#ifdef MAP_MIN_MEMORY
#include <map>
#endif

template<class Scalar> class FSCommPattern;

class CommunicatableObject
{
 protected:
  template<class Scalar>
    void setCommSize(FSCommPattern<Scalar> *pat, int localID, int remoteID);
  template<class Scalar>
    void pack(FSCommPattern<Scalar> *pat, int localID, int remoteID);
  template<class Scalar>
    void unpack(FSCommPattern<Scalar> *pat, int remoteID, int localID);
};

class GlobalToLocalMap : public CommunicatableObject
{
#ifdef MAP_MIN_MEMORY
  // HB: encapsulate a STL map<int,int>
  // I kept the min & max to optimized the operator[]; but maybe there already is
  // something similar in the find(...) method of the STL map.
  // A few optimizations can be done for the case where max-min+1 = size(): we just need
  // to send size()+3 data instead of 2.size()+3 data. If more, max-min+1 = size() AND min=0,
  // we may use a simple integer array to store globalToLocalArray as in this case the map
  // is merely permutation.
  private:
    int min, max; 
    std::map<int,int> globalToLocalArray;
#else
  private:
    int *globalToLocalArray;
    int min, max;
    int nKeys;
#endif  
  public:  
    GlobalToLocalMap();
    GlobalToLocalMap(int localSize, int *localToGlobalArray);
    ~GlobalToLocalMap();

    void initialize(int localSize, int *localToGlobalArray);
    /** \brief Obtain the localnumber or -1 if it was not found. */
    int operator[](int glNum) const;
    /** \brief Obtain the localnumber or throw if it was not found. */
    int at(int glNum) const { return globalToLocalArray.at(glNum); }

    void print(bool skip=false);
    int size();
    int numKeys(); 

    template<class Scalar>
      void setCommSize(FSCommPattern<Scalar> *pat, int localID, int remoteID);
    template<class Scalar>
      void pack(FSCommPattern<Scalar> *pat, int localID, int remoteID);
    template<class Scalar>
      void unpack(FSCommPattern<Scalar> *pat, int remoteID, int localID);
};

#endif
