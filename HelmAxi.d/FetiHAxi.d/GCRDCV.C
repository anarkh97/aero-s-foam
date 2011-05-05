#include <cmath>
#include <cstdio>
#include <Feti.d/DistrVector.h>
#include <HelmAxi.d/FetiHAxi.d/DistrComplexVector.h>
#include <HelmAxi.d/FetiHAxi.d/GCRDCV.h>


GCRDCV::GCRDCV(int maxSize, DistrInfo & dI) {

  int i;
  numP = 0;
  maxP = maxSize;

  // Copy construct distr info

  distrInfo = dI;
  distrInfo.domLen = new int[dI.numDom];
  for (i=0; i<dI.numDom; i++) distrInfo.domLen[i] = dI.domLen[i];

  // Initialize pointers 

  allFPFP = new DComplex[maxP];
  allP = new DistrComplexVector*[maxP];
  allFP = new DistrComplexVector*[maxP]; 

}


GCRDCV::~GCRDCV() {

  int i;

  for (i=0;i<numP; i++) { 
    delete allP[i];
    delete allFP[i];
  }

  delete[] allFP;
  delete[] allP;
  delete[] allFPFP;
  delete[] distrInfo.domLen;

}


void GCRDCV::orthoAdd(DistrComplexVector &p, DistrComplexVector &Fp, 
             DComplex FpFp) {

  // Vector assumed to be orthogonal

  if (numP == maxP) {
    numP--;
  } 
  else {
    allP[numP] = new DistrComplexVector(distrInfo);
    allFP[numP] = new DistrComplexVector(distrInfo);
  }

  // Copy the input into the set

  *allP[numP] = p;
  *allFP[numP] = Fp;
  allFPFP[numP] = FpFp;
  numP++;

}


void GCRDCV::orthogonalize(DistrComplexVector &r, DistrComplexVector &Fr, 
             DistrComplexVector &p, DistrComplexVector &Fp) {

  // r - the input vector to be orthogonalized to the current set
  // p - output vector 

  int i;
  p = r;
  Fp = Fr;

  if (numP==0) return;

  for(i=0; i<numP; i++) { 
    DComplex alpha = *allFP[i] ^ Fp;
    alpha = - alpha / allFPFP[i];
    p.linAdd(alpha, *allP[i]);
    Fp.linAdd(alpha, *allFP[i]);
  } 

}


int GCRDCV::predict(DistrComplexVector &r, DistrComplexVector &lambda0) {

  // r - the right hand side 
  // lambda0 - the modified initial solution

  if (numP == 0) return -1;

  int i;
  DComplex *y = new DComplex[numP];
  for (i=0; i<numP; i++) {
     y[i] = *allFP[i] ^ r;
     y[i] /= allFPFP[i];
     y[i] = - y[i];
  }

  for (i=0; i<numP; i++) 
      lambda0.linAdd(y[i], *allP[i]);

  delete[] y;
  return 0;

}


int GCRDCV::update(DistrComplexVector &x, DistrComplexVector &lambda0, 
            DistrComplexVector &r0) {

  if (numP == 0) return -1;

  int i;
  DComplex *y = new DComplex[numP];

  for (i=0; i<numP; i++) {
     y[i] = *allFP[i] ^ r0;
     y[i] /= allFPFP[i];
     y[i] = - y[i];
  }

  x = lambda0;

  for (i=0; i<numP; i++) 
       x.linAdd(y[i], *allP[i]);

  delete[] y;
  return 0;

}
