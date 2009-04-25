#ifndef _DISTRCOMPLEXVECTORSET_H_
#define _DISTRCOMPLEXVECTORSET_H_

#include <Utils.d/MyComplex.h>

class DistrComplexVector;
class DistrInfo;

class DistrComplexVectorSet {
  int    numVectors;
  DistrComplexVector *vecSet;

 public:

  // Constructor
  DistrComplexVectorSet(int numVec, DistrInfo &Info);

  DistrComplexVector& operator [] (int i) const;

  int numVec() { return numVectors; }
  double normTotal();
  DComplex ** globalSubData(int i);
  void zero();

};

#endif
