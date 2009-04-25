#ifndef _GLOBALTOLOCALMAP_H_
#define _GLOBALTOLOCALMAP_H_

//#define HB_USE_MEMCPY
#ifdef HB_USE_MEMCPY
#include <string.h>
#endif

//#define MAP_MIN_MEMORY
#ifdef MAP_MIN_MEMORY
#include <map>
#endif

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
  // PJSA: 2-14-02, this class is used for the globalToLocalMPC mapping
  // replaces previous integer array of size glMPCNum for every subdomain
/*  private:
    int localSize;
    int *localToGlobalArray; // must be sorted in ascending order;
    bool myData;
  public:
    GlobalToLocalMap() { localSize = 0; localToGlobalArray = 0; myData = false; }
    GlobalToLocalMap(int _localSize, int *_localToGlobalArray) {
      initialize(_localSize, _localToGlobalArray);
    }
    ~GlobalToLocalMap() { if(myData) delete [] localToGlobalArray; }
    void initialize(int _localSize, int *_localToGlobalArray) {
      localSize = _localSize;
      localToGlobalArray = _localToGlobalArray;
      myData = false;
    }
    int operator[](int glNum) { // binary search
      int low = 0, high = localSize-1;
      int middle;
      while(low <= high) {
        middle = (low+high)/2;
        if(localToGlobalArray[middle] == glNum) return middle;
        else if(localToGlobalArray[middle] < glNum) low = middle + 1;
        else high = middle - 1;
      }
      return -1;  // not found
    }
*/
//HB: encapsulate a STL map<int,int>
//I kept the min & max to optimized the operator[]; but maybe there already is
//something similar in the find(...) method of the STL map ...
//A few optimizations can be done for the case where max-min+1 = size(): we just need
//to send size()+3 data instead of 2.size()+3 data. If more, max-min+1 = size() AND min=0,
//we may use a simple integer array to store globalToLocalArray as in this case the map
//is merely permutation ...
  private:
    int min, max; 
    std::map<int,int> globalToLocalArray;
  
  public:  
    GlobalToLocalMap():globalToLocalArray() { min = max = 0; }
    GlobalToLocalMap(int localSize, int *localToGlobalArray) {
      initialize(localSize, localToGlobalArray);
    }
    ~GlobalToLocalMap() { if(~globalToLocalArray.empty()) globalToLocalArray.clear(); min = max = 0; }
    void initialize(int localSize, int *localToGlobalArray) {
      if(~globalToLocalArray.empty()) globalToLocalArray.clear();
      if(localSize<=0) { min = max = 0; return; }
      min = max = localToGlobalArray[0];
      if(min<0) { min = max = 0; }
      for(int i=0; i<localSize; ++i){
        if(localToGlobalArray[i]<0) continue; //HB: skip <0 values (for instance, -1 used as flag)
        if(localToGlobalArray[i] < min) min = localToGlobalArray[i];
        if(localToGlobalArray[i] > max) max = localToGlobalArray[i];
        globalToLocalArray[localToGlobalArray[i]] = i;
      }
    }

    int operator[](int glNum) {
      if((glNum < min) || (glNum > max)) return -1;
      else {
        std::map<int,int>::iterator I = globalToLocalArray.find(glNum);
        if(I!=globalToLocalArray.end()) return (*I).second;
        else return -1;
      }
    }

    void print(bool skip=false) {
      cerr << "GlobalToLocalMap of size "<<size()<<": ";
      std::map<int,int>::iterator I = globalToLocalArray.begin();
      while(I!=globalToLocalArray.end()){
        cerr <<"("<<(*I).first<<","<<(*I).second<<") ; ";
        I++;
      }
      cerr << endl;
    }

    int size() { return globalToLocalArray.size(); }
    int numKeys(){ return size(); } 
    template<class Scalar>
      void setCommSize(FSCommPattern<Scalar> *pat, int localID, int remoteID) {
        pat->setLen(localID, remoteID, 2*size()+3);
      }
    template<class Scalar>
      void pack(FSCommPattern<Scalar> *pat, int localID, int remoteID) {
        FSSubRecInfo<Scalar> sInfo = pat->getSendBuffer(localID, remoteID);
        int l = size();
        sInfo.data[0] = min;
        sInfo.data[1] = max;
        sInfo.data[2] = l;
        int i = 0;
        std::map<int,int>::iterator I = globalToLocalArray.begin();
        while(I!=globalToLocalArray.end()){
          sInfo.data[i+3]   = (*I).first; 
          sInfo.data[i+l+3] = (*I).second; 
          I++; i++;
        }
      }
    template<class Scalar>
      void unpack(FSCommPattern<Scalar> *pat, int remoteID, int localID) {
        if(~globalToLocalArray.empty()) globalToLocalArray.clear();
        FSSubRecInfo<Scalar> rInfo = pat->recData(remoteID, localID);
        min = rInfo.data[0];
        max = rInfo.data[1];
        int l = rInfo.data[2];
        for(int i=0; i<l; i++)
          globalToLocalArray[rInfo.data[i+3]] = rInfo.data[i+l+3];
      }
#else
  private:
    int *globalToLocalArray;
    int min, max;
    int nKeys;

  public:
    GlobalToLocalMap() { min = 0; max = -1; globalToLocalArray = 0; nKeys = 0; }
    GlobalToLocalMap(int localSize, int *localToGlobalArray) {
      initialize(localSize, localToGlobalArray);
    }
    ~GlobalToLocalMap() { if(globalToLocalArray) delete [] globalToLocalArray; }
    void initialize(int localSize, int *localToGlobalArray) {
      if(globalToLocalArray) delete [] globalToLocalArray;
      if(localSize<=0) { min = 0; max = -1; globalToLocalArray = 0; return; }
      min = max = localToGlobalArray[0];
      if(min<0) { min = max = 0; } //HB
      for(int i=0; i<localSize; ++i) {
        if(localToGlobalArray[i]<0) continue; //HB: skip <0 values (for instance, -1 used as flag)
        if(localToGlobalArray[i] < min) min = localToGlobalArray[i];
        if(localToGlobalArray[i] > max) max = localToGlobalArray[i];
      }
      if(min<0)
        cerr<<" ### PB in GlobalToLocalMap::initialize(...): min ("<<min<<") < 0"<<endl;
      globalToLocalArray = new int[max-min+1];
      for(int i=0; i < (max-min+1); ++i) globalToLocalArray[i] = -1;
      nKeys = 0;
      for(int i=0; i<localSize; ++i)
        if(localToGlobalArray[i]>=0) { 
          globalToLocalArray[localToGlobalArray[i]-min] = i;
          nKeys++;
        }
    }
    int operator[](int glNum) {
      if((glNum < min) || (glNum > max)) return -1;
      else return globalToLocalArray[glNum-min];
    }
    void print(bool skip=false) {
      cerr << "GlobalToLocalMap of size "<<max-min+1<<": ";
      for(int i=min; i<=max; ++i){
        if(skip && (*this)[i]<0) continue;
        cerr <<"("<<i<<","<<(*this)[i]<<") ; ";
      }
      cerr << endl;
    }
    int size() { return (max-min+1>0) ? max-min+1 : 0; } 
    int numKeys() { return nKeys; }
    long long mem() { return size()+3; }
    template<class Scalar>
      void setCommSize(FSCommPattern<Scalar> *pat, int localID, int remoteID) {
        pat->setLen(localID, remoteID, size()+3);
      }
    template<class Scalar>
      void pack(FSCommPattern<Scalar> *pat, int localID, int remoteID) {
        FSSubRecInfo<Scalar> sInfo = pat->getSendBuffer(localID, remoteID);
        sInfo.data[0] = min;
        sInfo.data[1] = max;
        sInfo.data[2] = nKeys;
#ifdef HB_USE_MEMCPY
        memcpy(&sInfo.data[3], globalToLocalArray, size()); //HB
#else
        for(int i=0; i<size(); ++i) sInfo.data[i+3] = globalToLocalArray[i];
#endif
      }
    template<class Scalar>
      void unpack(FSCommPattern<Scalar> *pat, int remoteID, int localID) {
        if(globalToLocalArray) delete [] globalToLocalArray;
        FSSubRecInfo<Scalar> rInfo = pat->recData(remoteID, localID);
        min  = rInfo.data[0];
        max  = rInfo.data[1];
        nKeys= rInfo.data[2];
        //cerr<<" in GlobalToLocalMap::unpack(...): min = "<<min<<", max = "<<max<<endl;
        if(max-min+1<=0) { min = 0; max = -1; globalToLocalArray = 0; nKeys = 0; return; }
        globalToLocalArray = new int[max-min+1];
#ifdef HB_USE_MEMCPY
        memcpy(globalToLocalArray,&rInfo.data[3], size()); //HB
#else
        for(int i=0; i<size(); ++i)  globalToLocalArray[i] = rInfo.data[i+3];
#endif
      }
#endif
};

#endif
