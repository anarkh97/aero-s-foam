#include <cstdio>

#include <Feti.d/DistrVector.h>
#include <HelmAxi.d/FetiHAxi.d/DistrComplexVectorSet.h>
#include <HelmAxi.d/FetiHAxi.d/DistrComplexVector.h>


DistrComplexVectorSet::DistrComplexVectorSet(int _numVectors, 
                       DistrInfo &Info)
{
 // This constructor creates a vector set of numVectors
 // each intialized to zero

 numVectors = _numVectors;
 DistrComplexVector *vecSet = new DistrComplexVector[numVectors];

 for (int i=0; i<numVectors; ++i) {
    DistrComplexVector v1(Info);
    vecSet[i].copy(v1);
 }

}


DistrComplexVector&
DistrComplexVectorSet::operator [] (int i) const
{
 return vecSet[i];
} 


double
DistrComplexVectorSet::normTotal() {

  int i;
  double norm=0.0;

  for (i=0; i<numVectors; ++i) 
     norm += real(vecSet[i]^vecSet[i]);     

  return norm; 

}


DComplex **
DistrComplexVectorSet::globalSubData(int i) {

  DComplex **buffer = new DComplex *[numVectors];
  for (int j=0; j<numVectors; ++j)  
     buffer[j] = vecSet[j].subData(i);

  return buffer;
}


void
DistrComplexVectorSet::zero() {

  for (int j=0; j<numVectors; ++j)
     vecSet[j].zero();

}


