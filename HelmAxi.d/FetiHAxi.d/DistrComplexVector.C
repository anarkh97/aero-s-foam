#include <Utils.d/dbg_alloca.h>
#include <math.h>
#include <Feti.d/DistrVector.h>
#include <Utils.d/MyComplex.h>
#include <HelmAxi.d/FetiHAxi.d/DistrComplexVector.h>
#include <Threads.d/Paral.h>


class VecOpComplex : public TaskDescr {
     DistrComplexVector *v1;
     DistrComplexVector *v2;
     DistrComplexVector *v3;
     DComplex *res;
     DComplex c;
     DComplex c1;
     DComplex c2;
     void (VecOpComplex::*f)(int);
   public:
     VecOpComplex(void (VecOpComplex::*_f)(int),
          DistrComplexVector *_v1, DistrComplexVector *_v2=0);

     VecOpComplex(void (VecOpComplex::*_f)(int), 
          DistrComplexVector *_v1, DistrComplexVector *_v2, DComplex c);

     VecOpComplex(void (VecOpComplex::*_f)(int),
          DistrComplexVector *_v1, DComplex c);

     VecOpComplex(void (VecOpComplex::*_f)(int), 
          DistrComplexVector *_v1, DistrComplexVector *_v2, DComplex *rc);

     VecOpComplex(void (VecOpComplex::*_f)(int),
          DistrComplexVector *_v1, DistrComplexVector *_v2, DComplex c,  
          DistrComplexVector *_v3);

     VecOpComplex(void (VecOpComplex::*_f)(int),
          DistrComplexVector *_v1, DComplex c1, DistrComplexVector *_v2, 
          DComplex c2,  DistrComplexVector *_v3);

     void alloc(int);
     void dotR(int);
     void dotC(int);
     void zero(int);
     void negate(int);
     void linAdd(int);
     void linAdd2(int);
     void linC(int);
     void linC2(int);
     void assign(int);
     void assign_times(int);
     void assign_plus(int);
     void assign_minus(int);
     void runFor(int);
};


VecOpComplex::VecOpComplex(void (VecOpComplex::*_f)(int), 
              DistrComplexVector *_v1, DistrComplexVector *_v2)
{
 f = _f;
 v1 = _v1;
 v2 = _v2;
}


VecOpComplex::VecOpComplex(void (VecOpComplex::*_f)(int), 
              DistrComplexVector *_v1, DistrComplexVector *_v2, DComplex _c)
{
 f = _f;
 v1 = _v1;
 v2 = _v2;
 c = _c;
}


VecOpComplex::VecOpComplex(void (VecOpComplex::*_f)(int), 
              DistrComplexVector *_v1, DComplex _c)
{
 f = _f;
 v1 = _v1;
 c = _c;
}


VecOpComplex::VecOpComplex(void (VecOpComplex::*_f)(int), 
              DistrComplexVector *_v1, DistrComplexVector *_v2,
              DComplex _c, DistrComplexVector *_v3)
{
 f = _f;
 v1 = _v1;
 v2 = _v2;
 c = _c;
 v3 = _v3;
}


VecOpComplex::VecOpComplex(void (VecOpComplex::*_f)(int), 
              DistrComplexVector *_v1, DComplex _c1, DistrComplexVector *_v2, 
              DComplex _c2, DistrComplexVector *_v3)
{
 f = _f;
 v1 = _v1;
 c1 = _c1;
 v2 = _v2;
 c2 = _c2;
 v3 = _v3;
}


VecOpComplex::VecOpComplex(void (VecOpComplex::*_f)(int), 
              DistrComplexVector *_v1, DistrComplexVector *_v2,
              DComplex *_r)
{
 f = _f;
 v1 = _v1;
 v2 = _v2;
 res = _r;
}


void
VecOpComplex::dotR(int threadNum)
{
 DComplex *d1 = v1->threadData(threadNum);
 DComplex *d2 = v2->threadData(threadNum);
 int len = v1->threadLen(threadNum);
 DComplex r = DComplex(0.0, 0.0);
 for(int i = 0; i < len; ++i)
    r += d1[i]*d2[i];
 res[threadNum] = r;
}


void
VecOpComplex::dotC(int threadNum)
{
 DComplex *d1 = v1->threadData(threadNum);
 DComplex *d2 = v2->threadData(threadNum);
 int len = v1->threadLen(threadNum);
 DComplex r = DComplex(0.0, 0.0);
 for(int i = 0; i < len; ++i)
    r += conj(d1[i])*d2[i];
 res[threadNum] = r;
}


void
VecOpComplex::linAdd(int threadNum)
{
 DComplex *d1 = v1->threadData(threadNum);
 DComplex *d2 = v2->threadData(threadNum);
 int len = v1->threadLen(threadNum);
 for(int i = 0; i < len; ++i)
    d1[i] += c*d2[i];
}


void
VecOpComplex::linAdd2(int threadNum)
{
 DComplex *d1 = v1->threadData(threadNum);
 DComplex *d2 = v2->threadData(threadNum);
 DComplex *d3 = v3->threadData(threadNum);
 int len = v1->threadLen(threadNum);
 for(int i = 0; i < len; ++i)
    d1[i] += c1*d2[i] + c2*d3[i];
}


void
VecOpComplex::linC(int threadNum)
{
 DComplex *d1 = v1->threadData(threadNum);
 DComplex *d2 = v2->threadData(threadNum);
 DComplex *d3 = v3->threadData(threadNum);
 int len = v1->threadLen(threadNum);
 for(int i = 0; i < len; ++i)
    d1[i] = d2[i] + c*d3[i];
}


void
VecOpComplex::linC2(int threadNum)
{
 DComplex *d1 = v1->threadData(threadNum);
 DComplex *d2 = v2->threadData(threadNum);
 DComplex *d3 = v3->threadData(threadNum);
 int len = v1->threadLen(threadNum);
 for(int i = 0; i < len; ++i)
    d1[i] = c1*d2[i] + c2*d3[i];
}


void
VecOpComplex::assign(int threadNum)
{
 DComplex *d1 = v1->threadData(threadNum);
 DComplex *d2 = v2->threadData(threadNum);
 int len = v1->threadLen(threadNum);
 for(int i = 0; i < len; ++i)
    d1[i] = d2[i];
}


void
VecOpComplex::assign_times(int threadNum)
{
 DComplex *d1 = v1->threadData(threadNum);
 int len = v1->threadLen(threadNum);
 for(int i = 0; i < len; ++i)
    d1[i] *= c;
}


void
VecOpComplex::assign_plus(int threadNum)
{
 DComplex *d1 = v1->threadData(threadNum);
 DComplex *d2 = v2->threadData(threadNum);
 int len = v1->threadLen(threadNum);
 for(int i = 0; i < len; ++i)
    d1[i] += d2[i];
}


void
VecOpComplex::assign_minus(int threadNum)
{
 DComplex *d1 = v1->threadData(threadNum);
 DComplex *d2 = v2->threadData(threadNum);
 int len = v1->threadLen(threadNum);
 for(int i = 0; i < len; ++i)
    d1[i] -= d2[i];
}


void
VecOpComplex::zero(int threadNum)
{
 DComplex *d1 = v1->threadData(threadNum);
 int len = v1->threadLen(threadNum);
 for(int i = 0; i < len; ++i)
    d1[i] = DComplex(0.0, 0.0);
}


void
VecOpComplex::negate(int threadNum)
{
 DComplex *d1 = v1->threadData(threadNum);
 int len = v1->threadLen(threadNum);
 int i;
 for(i=0; i<len; ++i)
   d1[i] = -d1[i];
}


void
VecOpComplex::alloc(int threadNum)
{
 v1->thV[threadNum] = new DComplex[v1->threadLen(threadNum)];
}


void
VecOpComplex::runFor(int threadNum)
{
 (this->*f)(threadNum);
}


DComplex
DistrComplexVector::ident()
{
 DComplex res = 0;
 for(int i=0; i < len; ++i)
   res += DComplex(i % 7)*v[i];
 return res;
}


DistrComplexVector::DistrComplexVector() {

 len = 0;
 numDom = 0;
 v = 0;
 subV = 0;
 subVLen = 0;
 nT = 0;
 thLen = 0;
 thV = 0;

}


DistrComplexVector::DistrComplexVector(const DistrInfo &info)
{
 // Set the overall length
 len = info.len;

 // Set number of Domains
 numDom = info.numDom;

 // Allocate memory for the values
 v = new DComplex[len];

 //array of double pointers to hold sub vector values
 subV = new DComplex *[info.numDom];

 // int array to hold sub vector lengths
 subVLen = new int[info.numDom];

 subV[0] = v;
 DComplex *v2 = v;

 nT = threadManager->numThr();
 int iThread, md;

 thV   = new DComplex *[nT];
 thLen = new int[nT];

 int tLen;
 for(iThread = 0; iThread < nT; ++iThread) {
   tLen = 0;
   thV[iThread] = v2;
   for(md = iThread; md < info.numDom; md += nT) {
     subV[md] = v2;
     subVLen[md] = info.domLen[md];
     v2 += info.domLen[md];
     tLen += info.domLen[md];
   }
   thLen[iThread] = tLen;
 }

 /*
 VecOp makeAlloc(VecOp::alloc, this);
 threadManager->execParal(nT, &makeAlloc);
 for(iThread = 0; iThread < nT; ++iThread) {
   DComplex *v2 = thV[iThread];
   for(md = iThread; md < info.numDom; md += nT) {
      subV[md] = v2;
      v2 += info.domLen[md];
   }
 }
 */

}


void
DistrComplexVector::zero()
{
 VecOpComplex zeroAll(&VecOpComplex::zero,this);
 threadManager->execParal(nT, &zeroAll);
}


DistrComplexVector &
DistrComplexVector::copy( DistrComplexVector &x) {

 len     = x.len;
 numDom  = x.numDom;
 v       = x.v;
 subV    = x.subV;
 subVLen = x.subVLen;
 nT      = x.nT;
 thLen   = x.thLen;
 thV     = x.thV;

 return *this;

}


void
DistrComplexVector::negate()
{
 VecOpComplex negateAll(&VecOpComplex::negate,this);
 threadManager->execParal(nT, &negateAll);
}


DComplex 
DistrComplexVector::norm()
{
 return sqrt((*this)^(*this));
}


DComplex
DistrComplexVector::operator * (DistrComplexVector&x)
{

 DComplex *partial = (DComplex *) dbg_alloca(sizeof(DComplex)*nT);
 VecOpComplex dotAll(&VecOpComplex::dotR, this, &x, partial);
 threadManager->execParal(nT,&dotAll);
 DComplex res = DComplex(0.0, 0.0);
 for(int i=0; i < nT; ++i)
    res += partial[i];

 return res;
}


DComplex
DistrComplexVector::operator ^ (DistrComplexVector&x)
{

 DComplex *partial = (DComplex *) dbg_alloca(sizeof(DComplex)*nT);
 VecOpComplex dotAll(&VecOpComplex::dotC, this, &x, partial);
 threadManager->execParal(nT,&dotAll);
 DComplex res = DComplex(0.0, 0.0);
 for(int i=0; i < nT; ++i) {
    res += partial[i];
 }

 return res;
}


#include <stdio.h>

DistrComplexVector &
DistrComplexVector::operator=(DistrComplexVector &x)
{
 if(x.len != len) {
  fprintf(stderr, "Length error in = %d not equal to %d\n",x.len, len);
 }
 VecOpComplex assign(&VecOpComplex::assign,this, &x);
 threadManager->execParal(nT, &assign);
 return *this;
}


DistrComplexVector &
DistrComplexVector::operator*=(DComplex c)
{
 VecOpComplex assign_times(&VecOpComplex::assign_times,this, c);
 threadManager->execParal(nT, &assign_times);
 return *this;
}


DistrComplexVector &
DistrComplexVector::operator+=(DistrComplexVector &x)
{
 if(x.len != len) {
  fprintf(stderr, "Length error in +=\n");
 }
 VecOpComplex assign_plus(&VecOpComplex::assign_plus,this, &x);
 threadManager->execParal(nT, &assign_plus);
 return *this;
}


DistrComplexVector &
DistrComplexVector::operator-=(DistrComplexVector &x)
{
 if(x.len != len) {
  fprintf(stderr, "Length error in -=\n");
 }
 VecOpComplex assign_minus(&VecOpComplex::assign_minus,this, &x);
 threadManager->execParal(nT, &assign_minus);
 return *this;
}


DistrComplexVector &
DistrComplexVector::linAdd(DComplex c, DistrComplexVector &x)
{
 if(x.len != len) {
  fprintf(stderr, "Length error in linAdd\n");
 }
 VecOpComplex add(&VecOpComplex::linAdd, this, &x, c);
 threadManager->execParal(nT, &add);
 return *this;
}


DistrComplexVector &
DistrComplexVector::linAdd(DComplex c1, DistrComplexVector&x, DComplex c2, DistrComplexVector&y)
{
 if(x.len != y.len) {
  fprintf(stderr, "Length error in linAdd\n");
 }
 VecOpComplex add2(&VecOpComplex::linAdd2, this, c1, &x, c2, &y);
 threadManager->execParal(nT, &add2);
 return *this;
}


DistrComplexVector &
DistrComplexVector::linC( DistrComplexVector &x, DComplex c, DistrComplexVector &y)
{
 if(x.len != y.len) {
  fprintf(stderr, "Length error in linC\n");
 }
 VecOpComplex addC(&VecOpComplex::linC, this, &x, c, &y);
 threadManager->execParal(nT, &addC);
 return *this;
}


DistrComplexVector &
DistrComplexVector::linC( DComplex c1, DistrComplexVector &x, DComplex c2, DistrComplexVector &y)
{
 if(x.len != y.len) {
  fprintf(stderr, "Length error in linC\n");
 }
 VecOpComplex linC2(&VecOpComplex::linC2, this, c1, &x, c2, &y);
 threadManager->execParal(nT, &linC2);
 return *this;
}


#include <stdio.h>
void
DistrComplexVector::print() {

 for(int i=0; i < len; ++i)
   fprintf(stderr," v[%d] = (%e,%e) \n",i,real(v[i]), imag(v[i]));
 fprintf(stderr,"\n");

}


void
DistrComplexVector::printAll() {

 fprintf(stderr,"Length of DistrComplexVector = %d\n",len);
 fprintf(stderr,"Number of Domains = %d\n",numDom);
 int i,j;
 for(i=0; i<numDom; ++i) {
   fprintf(stderr,"--- Sub Vector Length = %d\n",subVLen[i]);
   for(j=0; j<subVLen[i]; ++j)
     fprintf(stderr,"v(%d) = (%e,%e) \n",j+1,real(subV[i][j]),
                    imag(subV[i][j]));
 }

}

