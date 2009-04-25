#include <math.h>
#include <stdio.h>
#include <Utils.d/dbg_alloca.h>
#include <Utils.d/linkfc.h>
#include <Threads.d/Paral.h>
#include <Timers.d/DistTimer.h>
#include <Timers.d/GetTime.h>
#include <HelmAxi.d/FetiHAxi.d/GCRC.h>


extern "C" {
void _FORTRAN(zgemv)(const char &, const int &, const int &, const DComplex &,
                 DComplex *, const int &, DComplex *, const int &,
                 const DComplex &, DComplex *, const int &);
};


GCRC::GCRC(int _len, int maxSize) {

 len = _len;
 maxP = maxSize;
 numP = 0;
 allFPFiP = new DComplex[maxP];

 numTasks = threadManager->numThr();
 oos = new TaskDescr*[numTasks];
 int locLen, index=0;
 int avl = len/numTasks;
 int mod = len % numTasks;

 for(int i = 0; i < numTasks; ++i) {
    locLen = (i < mod) ? avl+1 : avl;
    oos[i] = new GCRCOp(this, index, locLen);
    index += locLen;
 }

}


void
GCRC::orthoAdd(DComplex *p, DComplex *Fp, DComplex FpFp) {

 if (numP == maxP) 
     numP--;

 allFPFiP[numP] = FpFp;
 numP ++;

 // op1 = p.data();
 // op2 = Fp.data();
 op1 = p;
 op2 = Fp;
 operation = &GCRCOp::addVec;
 threadManager->execParal(numTasks, oos);

}


void
GCRC::orthoAddTimed(DistTimer &timer, DComplex *p, DComplex *Fp, 
                    DComplex FpFp) {

 double initTime = getTime();
 long initMem  = threadManager->memoryUsed();

 if (numP == maxP) 
     numP--;

 allFPFiP[numP] = FpFp;
 numP ++;

 // op1 = p.data();
 // op2 = Fp.data();
 op1 = p;
 op2 = Fp;
 operation = &GCRCOp::addVec;
 threadManager->execTimedParal(timer, numTasks, oos);

 timer.addOverAll( threadManager->memoryUsed()-initMem, getTime()-initTime );

}


void
GCRC::orthogonalize(DComplex *r, DComplex *Fr, DComplex *p, DComplex *Fp) {

  int i;
  DComplex *y = (DComplex *)dbg_alloca(numP*sizeof(DComplex));
   
  op2 = y;
  op1 = Fr;
  op3 = p;
  for(i = 0; i < numP; ++i)
     y[i] = DComplex(0.0, 0.0);
  operation = &GCRCOp::Fdot;
  threadManager->execParal(numTasks, oos);

  for(i = 0; i < numP; ++i) 
     y[i] /= allFPFiP[i];

  op1 = r;
  op3 = p;
  operation = &GCRCOp::multAdd;
  threadManager->execParal(numTasks, oos);

  op1 = Fr;
  op3 = Fp;
  operation = &GCRCOp::multFAdd;
  threadManager->execParal(numTasks, oos);

}


void
GCRC::orthogonalizeTimed(DistTimer &timer, DComplex *r, DComplex *Fr, 
                         DComplex *p, DComplex *Fp) {

  double initTime = getTime();
  long initMem = threadManager->memoryUsed();

  int i;
  DComplex *y = (DComplex *)dbg_alloca(numP*sizeof(DComplex));
   
  op2 = y;
  op1 = Fr;
  op3 = p;
  for(i = 0; i < numP; ++i)
     y[i] = DComplex(0.0, 0.0);
  operation = &GCRCOp::Fdot;
  threadManager->execTimedParal(timer, numTasks, oos);

  for(i = 0; i < numP; ++i) 
     y[i] /= allFPFiP[i];

  op1 = r;
  op3 = p;
  operation = &GCRCOp::multAdd;
  threadManager->execTimedParal(timer, numTasks, oos);

  op1 = Fr;
  op3 = Fp;
  operation = &GCRCOp::multFAdd;
  threadManager->execTimedParal(timer, numTasks, oos);

  timer.addOverAll( threadManager->memoryUsed()-initMem, getTime()-initTime);

}


void
GCRC::predict(DComplex *r, DComplex *lambda) {

  int i;
  DComplex *y = (DComplex *)dbg_alloca(numP*sizeof(DComplex));


  op2 = y;
  op1 = r; // .data();
  op3 = lambda; // .data();
  for(i = 0; i < numP; ++i)
      y[i] = 0;
  operation = &GCRCOp::Fdot;
  threadManager->execParal(numTasks, oos);

  for(i = 0; i < numP; ++i)
     y[i] /= allFPFiP[i];

  operation = &GCRCOp::mult;
  threadManager->execParal(numTasks, oos);

}

GCRCOp::GCRCOp(GCRC *_os, int _index, int _loclen) {

 os = _os; index = _index; loclen = _loclen;
 locAllP = locAllFiP = 0;
 numP = 0;

}

void
GCRCOp::run() {

 (this->*os->operation)();

}

void
GCRCOp::addVec() {

 if(numP == os->maxP)
   numP--;

 if(locAllP == 0) {
   locAllP   = new DComplex[2*os->maxP*loclen];
   locAllFiP = locAllP + os->maxP*loclen;
 }

 for(int i = 0; i < loclen; ++ i) {
   locAllP[i+ numP*loclen]   = os->op1[i+index];
   locAllFiP[i+ numP*loclen] = os->op2[i+index];
 }

 numP++;

}


void
GCRCOp::Fdot() {

  int i;
  char trans = 'C';
  DComplex *y = (DComplex *)dbg_alloca(os->maxP*sizeof(DComplex));

  _FORTRAN(zgemv)(trans, loclen, numP, DComplex(-1.0,0.0), locAllFiP,
         loclen, os->op1+index, 1, DComplex(0.0,0.0), y, 1);

  os->lock.lock();

  for(i = 0; i < os->numP; ++i)
     os->op2[i] += y[i];

  os->lock.unlock();

}

void
GCRCOp::dot() {

  int i;
  char trans = 'C';
  DComplex *y = (DComplex *)dbg_alloca(os->maxP*sizeof(DComplex));

  _FORTRAN(zgemv)(trans, loclen, numP, DComplex(-1.0,0.0), locAllP,
         loclen, os->op1+index, 1, DComplex(0.0,0.0), y, 1);

  os->lock.lock();

  for(i = 0; i < os->numP; ++i)
     os->op2[i] += y[i];

  os->lock.unlock();
}


void
GCRCOp::mult() {

  char trans = 'N';
 
  _FORTRAN(zgemv)(trans, loclen, numP, DComplex(1.0,0.0), locAllP,
         loclen,  os->op2, 1, DComplex(1.0,0.0), os->op3+index, 1);

}


void
GCRCOp::multAdd() {

  char trans = 'N';

  for(int i=0; i< loclen; ++i)
     os->op3[index+i] = os->op1[index+i];

  _FORTRAN(zgemv)(trans, loclen, numP, DComplex(1.0,0.0), locAllP,
         loclen,  os->op2, 1, DComplex(1.0,0.0), os->op3+index, 1);

}


void
GCRCOp::multFAdd() {

  char trans = 'N';

  for(int i=0; i< loclen; ++i)
     os->op3[index+i] = os->op1[index+i];

  _FORTRAN(zgemv)(trans, loclen, numP, DComplex(1.0,0.0), locAllFiP,
         loclen,  os->op2, 1, DComplex(1.0,0.0), os->op3+index, 1);

}
