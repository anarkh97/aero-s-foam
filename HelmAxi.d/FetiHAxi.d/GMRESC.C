#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <Utils.d/dbg_alloca.h>
#include <Utils.d/linkfc.h>
#include <Threads.d/Paral.h>
#include <HelmAxi.d/FetiHAxi.d/GMRESC.h>


#define H(i,j) matrixH[(j)*(maxV+1) + (i)]
#define v(i) (*allV[i])
#define w(i) (*allW[i])


extern "C" {
void _FORTRAN(zgemv)(const char &, const int &, const int &, 
                 const DComplex &, DComplex *, const int &, DComplex *,
                 const int &, const DComplex &, DComplex *, const int &);
};


void GMRESC::generateRotation(DComplex a, DComplex b, DComplex &cs, 
             DComplex &ss) {

  double temp = sqrt(norm(a) + norm(b));
  cs = a / temp;
  ss = conj(b) / temp; 

}


void GMRESC::applyRotation(DComplex &a, DComplex &b, DComplex cs, 
             DComplex ss) {

  DComplex temp1 = conj(cs) * a + ss * b;
  DComplex temp2 = -conj(ss) * a + cs * b;
  a = temp1;
  b = temp2;

}


GMRESC::GMRESC(int _len, int maxSize) {

  len = _len;
  numV = 0;
  maxV = maxSize;

  numTasks = threadManager->numThr();
  oos = new TaskDescr*[numTasks];
  int locLen, index=0;
  int avl = len/numTasks;
  int mod = len % numTasks;

  for(int i = 0; i < numTasks; ++i) {
    locLen = (i < mod) ? avl+1 : avl;
    oos[i] = new GMRESCOp(this, index, locLen);
    index += locLen;
  }

  DComplex * p = new DComplex[4*(maxV+1)];
  givensC = p;
  givensS = p + maxV+1;
  g = p + 2* (maxV+1);
  y = p + 3* (maxV+1);
  matrixH = new DComplex[(maxV+1)*maxV];
  reInit();

}


GMRESC::~GMRESC() {

  delete[] givensC;
  delete[] matrixH;

}


void GMRESC::reInit() {

  numV = 0;

  for(int i=0; i < maxV+1 ; i++) { 
    g[i] = DComplex(0.0,0.0); 
    givensC[i] = DComplex(0.0,0.0); 
    givensS[i] = DComplex(0.0,0.0); 
    for(int j=0; j < maxV ; j++) H(i,j) = DComplex(0.0,0.0);
  }

  operation = &GMRESCOp::reInit;
  threadManager->execParal(numTasks, oos);

}


void GMRESC::init(DComplex *v0, double beta) {

  g[0] = DComplex(beta,0);
  op1 = v0;
  operation = &GMRESCOp::addVec;
  threadManager->execParal(numTasks, oos);
  numV++;

}


double GMRESC::orthoAdd(DComplex *Fv, DComplex *v) {

  int i;

  if (numV == maxV+1) { 
    fprintf(stderr, "Maximum number of directions exceeeded in GMRES.\n");
    fprintf(stderr, "Terminating.\n");
    exit(-1);
  }

  op1 = Fv;
  op2 = y;

  for(i = 0; i < numV; ++i)
      y[i] = 0;

  operation = &GMRESCOp::dot;
  threadManager->execParal(numTasks, oos);

  for(i = 0; i < numV; ++i) {
     y[i] = DComplex(2.0,0.0) * y[i];
     H(i,numV-1) = y[i];
  }

  op3 = v;
  operation = &GMRESCOp::multAdd;
  threadManager->execParal(numTasks, oos);

  double nrm=0.0;
  for(i = 0; i < len; ++i) nrm += real(conj(v[i])*v[i]);
  H(numV,numV-1) = DComplex(sqrt(2.0*nrm),0.0);

  op2 = & H(numV,numV-1);
  op1 = v;
  operation = &GMRESCOp::addVecAndNorm;
  threadManager->execParal(numTasks, oos);

  for (i = 0; i < numV-1; i++)
     applyRotation(H(i,numV-1), H(i+1,numV-1), givensC[i], givensS[i]);

  generateRotation(H(numV-1,numV-1), H(numV,numV-1), givensC[numV-1], 
                   givensS[numV-1]);
  applyRotation(H(numV-1,numV-1), H(numV,numV-1), givensC[numV-1], 
                givensS[numV-1]);
  applyRotation(g[numV-1], g[numV], givensC[numV-1], givensS[numV-1]);

  numV++;
  return norm(g[numV-1]);

}


void GMRESC::solution(DComplex *u) {

  int i,j;

//Back substitute

  for(i=numV-2;i>=0;i--) {
    y[i] = g[i];
    for(j=numV-2; j>i;j--) y[i] -= H(i,j) * y[j];
    y[i] /= H(i,i);
  }

  op1 = u;
  op2 = y;
  operation = &GMRESCOp::mult;
  threadManager->execParal(numTasks, oos);

}


GMRESCOp::GMRESCOp(GMRESC *_os, int _index, int _loclen) {

 os = _os; index = _index; loclen = _loclen;
 locAllV = 0;
 numV = 0;

}


void GMRESCOp::reInit() {

 numV = 0;

}


void
GMRESCOp::run() {

 (this->*os->operation)();

}


void
GMRESCOp::addVec() {

 if(locAllV == 0) {
   locAllV   = new DComplex[(os->maxV+1)*loclen];
 }

 for(int i = 0; i < loclen; ++ i) {
   locAllV[i+ numV*loclen]   = os->op1[i+index];
 }

 numV++;

}


void
GMRESCOp::addVecAndNorm() {

 int i; 
 DComplex nrm = os->op2[0];

 for(i = 0; i < loclen; ++ i) {
   locAllV[i+ numV*loclen]   = os->op1[i+index] / nrm;
   os->op1[i+index] = locAllV[i+ numV*loclen];
 }

 numV++;

}


void
GMRESCOp::test() {

  int i,j;
  int x=0;

  fprintf(stderr,"numV:%d\n",numV);

  for(i = 0; i < numV; ++i) for(j=0;j<numV;j++) {
    DComplex inner(0.0,0.0);
    for(int k = 0; k < loclen; ++ k) {
       inner += conj(locAllV[k+ i*loclen]) * locAllV[k+ j*loclen];
    }
    os->lock.lock();
    os->op1[x] += inner;
    os->lock.unlock();
    x++;
  }

}


void
GMRESCOp::dot() {

  int i;
  char trans = 'C';
  DComplex *y = (DComplex *)dbg_alloca(os->maxV*sizeof(DComplex));

  _FORTRAN(zgemv)(trans, loclen, numV, DComplex(1.0,0.0), locAllV,
         loclen, os->op1+index, 1, DComplex(0.0,0.0), y, 1);
  os->lock.lock();

  for(i = 0; i < os->numV; ++i)
     os->op2[i] += y[i];

  os->lock.unlock();

}


void
GMRESCOp::mult() {

  char trans = 'N';
 
  _FORTRAN(zgemv)(trans, loclen, numV-1, DComplex(-1.0,0.0), locAllV,
         loclen,  os->op2, 1, DComplex(0.0,0.0), os->op1+index, 1);

}


void
GMRESCOp::multAdd() {

  char trans = 'N';

  for(int i=0; i< loclen; ++i)
     os->op3[index+i] = os->op1[index+i];

  _FORTRAN(zgemv)(trans, loclen, numV, DComplex(-1.0,0.0), locAllV,
         loclen,  os->op2, 1, DComplex(1.0,0.0), os->op3+index, 1);

}
