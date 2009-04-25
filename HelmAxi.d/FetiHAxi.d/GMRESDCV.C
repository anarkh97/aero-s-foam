#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <Feti.d/DistrVector.h>
#include <HelmAxi.d/FetiHAxi.d/DistrComplexVector.h>
#include <HelmAxi.d/FetiHAxi.d/GMRESDCV.h>


#define H(i,j) matrixH[(j)*(maxV+1) + i]
#define v(i) (*allV[i])
#define w(i) (*allW[i])


void GMRESDCV::generateRotation(DComplex a, DComplex b, DComplex &cs, 
               DComplex &ss) {

  double temp = sqrt(norm(a) + norm(b));
  cs = a / temp;
  ss = conj(b) / temp; 

}


void GMRESDCV::applyRotation(DComplex &a, DComplex &b, DComplex cs, 
               DComplex ss) {

  DComplex temp1 = conj(cs) * a + ss * b;
  DComplex temp2 = -conj(ss) * a + cs * b;
  a = temp1;
  b = temp2;

}


GMRESDCV::GMRESDCV(int maxSize, DistrInfo & dI) {

  int i,j;
  numV = 0;
  maxV = maxSize;

  // Copy construct distr info

  distrInfo = dI;
  distrInfo.domLen = new int[dI.numDom];
  for (i=0; i<dI.numDom; i++) 
       distrInfo.domLen[i] = dI.domLen[i];

  // Initialize pointers 

  allV = new DistrComplexVector*[maxV+1];
  allW = new DistrComplexVector*[maxV]; 
  givensC = new DComplex[maxV+1];
  givensS = new DComplex[maxV+1];
  g = new DComplex[maxV+1]; 
  y = new DComplex[maxV+1]; 
  matrixH = new DComplex[(maxV+1)*maxV];

  for(i=0; i < maxV+1 ; i++) { 
    g[i] = DComplex(0.0,0.0); 
    for(j=0; j < maxV ; j++) H(i,j) = DComplex(0.0,0.0);
  }

}


GMRESDCV::~GMRESDCV() {

  int i;

  for (i=0;i<numV; i++) { 
    delete allV[i];
    delete allW[i];
  }

  delete allV[numV];
  delete[] y;
  delete[] g;
  delete[] givensS;
  delete[] givensC;
  delete[] allW;
  delete[] allV;
  delete[] matrixH;
  delete[] distrInfo.domLen;

}


void GMRESDCV::reInit() {

  int i,j;
  for (i=0;i<numV; i++) { 
    delete allV[i];
    delete allW[i];
  }

  delete allV[numV];

  numV = 0;

  for (i=0; i < maxV+1 ; i++) { 
    g[i] = DComplex(0.0,0.0); 
    for (j=0; j < maxV ; j++) 
        H(i,j) = DComplex(0.0,0.0);
  }

}


void GMRESDCV::init(DistrComplexVector & r0) {

  beta = sqrt(real(r0^r0));
  allV[0] = new DistrComplexVector(distrInfo);
  v(0).zero();
  v(0).linAdd(DComplex(1.0/beta,0.0),r0);
  g[0] = DComplex(beta,0);
  r0 = v(0);

}


void GMRESDCV::orthoAdd(DistrComplexVector & Fv, DistrComplexVector & v) {

  int i;

  if (numV == maxV) { 
    fprintf(stderr, "Maximum number of directions exceeeded in GMRES.\n");
    fprintf(stderr, "Terminating.\n");
    exit(-1);
  } 

  allW[numV] = new DistrComplexVector(distrInfo);
  w(numV) = Fv;

  for(i=0; i<=numV; i++) {
    H(i,numV) =  v(i) ^ w(numV) ;
    w(numV).linAdd( -H(i,numV), v(i));
  }

  H(numV+1,numV) = DComplex(sqrt(real(w(numV)^w(numV))),0.0);
  allV[numV+1] = new DistrComplexVector(distrInfo);
  v(numV+1).zero();
  v(numV+1).linAdd(DComplex(1.0,0.0)/H(numV+1,numV),w(numV));

  v = v(numV+1);
//  dumpH();

}


double GMRESDCV::update() {

  int k;

  for (k = 0; k < numV; k++)
     applyRotation(H(k,numV), H(k+1,numV), givensC[k], givensS[k]);

  generateRotation(H(numV,numV), H(numV+1,numV), givensC[numV], givensS[numV]);
  applyRotation(H(numV,numV), H(numV+1,numV), givensC[numV], givensS[numV]);
  applyRotation(g[numV], g[numV+1], givensC[numV], givensS[numV]);

//  dumpH();
  numV++;

  return norm(g[numV]);

}


void GMRESDCV::solution(DistrComplexVector &u) {

  int i,j;

//Back substitute

  for(i=numV-1;i>=0;i--) {
    y[i] = g[i];
    for(j=numV-1; j>i;j--) y[i] -= H(i,j) * y[j];
    y[i] /= H(i,i);
  }

  u.zero();

  for(i=0;i<numV;i++) 
    u.linAdd(-y[i],v(i));


}


void GMRESDCV::dumpH() {

  int i,j;

  fprintf(stderr,"\n");

  for(i=0; i <= numV+1 ; i++) { 
    fprintf(stderr,"%.4e&%.4e  ",real(g[i]),imag(g[i]));
    for (j=0; j <= numV; j++) 
        fprintf(stderr,"%.4e&%.4e ",real(H(i,j)),imag(H(i,j)));
    fprintf(stderr,"\n");
  }

  fflush(stderr);

}
