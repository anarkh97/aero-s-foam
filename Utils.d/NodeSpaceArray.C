#include <Utils.d/NodeSpaceArray.h>
#include <stdio.h>
#include <stdlib.h>



Contraction
Tensor::operator | (Tensor &b) { return Contraction(*this, b); }

DoubleContraction 
Tensor::operator || (Tensor &b) { return DoubleContraction(*this, b); }


void 
Tensor::dblContractInto(const Tensor &, Tensor *) const {

fprintf(stderr," Be careful my dear, it seems that the operands don`t fit the operator ||...   \n");
     exit(-1);
}


void 
Tensor::splContractInto(const Tensor &, Tensor *) const {

fprintf(stderr," Be careful my dear, it seems that the operands don`t fit the operator |...   \n");
     exit(-1);
}




Tensor_d0s2
Tensor_d0s2::operator + (const Tensor_d0s2 &tens) const{
Tensor_d0s2 t;
  for (int i = 0; i < 9; i++)
    t[i]=v[i]+tens[i];
return t;
}

Tensor_d0s2
Tensor_d0s2::operator * (double scal){
Tensor_d0s2 t;
  for (int i = 0; i < 9; i++)
    t[i]=v[i]*scal;
return t;
}



Tensor_d0s2
operator *(double scal, const Tensor_d0s2 &tens){
Tensor_d0s2 t;
  for (int i = 0; i < 9; i++)
    t[i]=tens[i]*scal;
return t;
}







/*
Tensor_d0s2
Tensor_d0s2::operator | (const Tensor_d0s2 &tens) const {
Tensor_d0s2 t;
for (int i = 0; i < 3; i++)
  for (int j = 0; j < 3; j++)
    for (int k = 0; k < 3; k++)
      t[3*i+j]+=v[3*i+k]*tens[3*k+j];
return t;
}


Tensor_d1s2_full
Tensor_d0s2::operator | (const Tensor_d1s2_full &tens)const {
int size = tens.getSize();
Tensor_d1s2_full t(size);
for (int m = 0; m < size; m++)
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++)
      for (int k = 0; k < 3; k++)
        t[m][3*i+j]+=v[3*i+k]*tens[m][3*k+j];
return t;
}


Tensor_d1s2_full
Tensor_d0s2::operator | (const Tensor_d1s2_sparse &tens)const {
int size = tens.getSize();
Tensor_d1s2_full t(size);
for (int n = 0; n < size; n++)
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++)
      t[n][3*i+j]=v[3*i+(n%3)]*tens[(n-(n%3))/3][3*(n%3)+j];
return t;

}
*/


void
Tensor_d0s2::splContractInto(const Tensor &b, Tensor *result) const
{
 
const Tensor_d0s2 *tens = dynamic_cast<const Tensor_d0s2 *>(&b);
  if(tens)
  {
    Tensor_d0s2 &t = static_cast<Tensor_d0s2 &>(*result);  
    for (int i = 0; i < 3; i++)
      for (int j = 0; j < 3; j++){
        t[3*i+j]=0.0;
        for (int k = 0; k < 3; k++)
          t[3*i+j]+=v[3*i+k]*(*tens)[3*k+j];
                                  }
  }
  else
  {
    const Tensor_d1s2_sparse *tens1 = dynamic_cast<const Tensor_d1s2_sparse *>(&b);
    if(tens1)
    {
      Tensor_d1s2_full &t = static_cast<Tensor_d1s2_full &>(*result); 
      int size = tens1->getSize();
      for (int n = 0; n < size; n++)
        for (int i = 0; i < 3; i++)
          for (int j = 0; j < 3; j++)
            t[n][3*i+j]=v[3*i+(n%3)]*(*tens1)[(n-(n%3))/3][3*(n%3)+j];
    }
    else
    {
      const Tensor_d1s2_full &tens2 = static_cast<const Tensor_d1s2_full &>(b);
      Tensor_d1s2_full &t = static_cast<Tensor_d1s2_full &>(*result);
      int size = tens2.getSize();
      for (int m = 0; m < size; m++)
        for (int i = 0; i < 3; i++)
          for (int j = 0; j < 3; j++){
            t[m][3*i+j]=0.0;
            for (int k = 0; k < 3; k++)
              t[m][3*i+j]+=v[3*i+k]*tens2[m][3*k+j];
                                      }
    }
  }
}

void
Tensor_d0s2::getDeterminant(double &det)  {               det=v[0]*v[4]*v[8]+v[3]*v[7]*v[2]+v[1]*v[5]*v[6]-(v[6]*v[4]*v[2]+v[3]*v[1]*v[8]+v[0]*v[7]*v[5]);
}


void
Tensor_d0s2::getInverse(Tensor_d0s2 &t)  {  
double d;
  getDeterminant(d);
  if ( d==0.0 )
    {fprintf(stderr," Caution : This Matrix is singular  \n");
     exit(-1);
    }
  t[0]=(1/d)*(v[4]*v[8]-v[5]*v[7]);
  t[1]=(1/d)*(-v[1]*v[8]+v[2]*v[7]);
  t[2]=(1/d)*(v[1]*v[5]-v[2]*v[4]);
  
  t[3]=(1/d)*(v[6]*v[5]-v[3]*v[8]);
  t[4]=(1/d)*(-v[6]*v[2]+v[0]*v[8]);
  t[5]=(1/d)*(v[3]*v[2]-v[0]*v[5]);
  
  t[6]=(1/d)*(-v[6]*v[4]+v[3]*v[7]);
  t[7]=(1/d)*(v[6]*v[1]-v[0]*v[7]);
  t[8]=(1/d)*(-v[3]*v[1]+v[0]*v[4]);
}


void 
Tensor_d0s2::getTranspose(Tensor_d0s2 &t){
for (int i = 0; i < 3; i++)
  for (int j = 0; j < 3; j++)
    t[3*i+j]=v[3*j+i];
}


void 
Tensor_d0s2::convertToSym(Tensor_d0s2_Ss12 &t){
for (int i = 0; i < 3; i++)
  for (int j = i; j < 3; j++)
    t[i*(5-i)/2+j]=v[3*i+j];
}     














Tensor_d1s0::Tensor_d1s0(int _size){
size=_size;
v = new double[size];
for (int i = 0; i<size; i++) 
  v[i]=0;
}


Tensor_d1s0
Tensor_d1s0::operator + (const Tensor_d1s0 &tens)const {
Tensor_d1s0 t(size);
for (int k = 0; k < size; k++)
  t[k]=v[k]+tens[k];
return t;
}


Tensor_d1s0 
operator *(double d, const Tensor_d1s0 &tens){
Tensor_d1s0 t(tens.size);
for (int k = 0; k < t.size; k++)
  t[k]=d*tens[k];
return t;
};


Tensor_d1s0 &
Tensor_d1s0::operator=(const Tensor_d1s0 &t) {
if(size != t.size) {
  if(v) delete[] v;
  size = t.size;
  v = new double[size];
                    }
for(int i = 0; i < size; ++i)
  v[i] = t.v[i];
return *this;
}


Tensor_d1s0::Tensor_d1s0(const Tensor_d1s0 &t) {
size = t.size;
v = new double[size];
for(int i = 0; i < size; ++i)
  v[i] = t.v[i];
}
















Tensor_d1s1::Tensor_d1s1(int _size){
size=_size;
v = new Tensor_d0s1[size];
}


Tensor_d1s1 &
Tensor_d1s1::operator=(const Tensor_d1s1 &t) {
if(size != t.size) {
  if(v) delete[] v;
  size = t.size;
  v = new Tensor_d0s1[size];
                   }
for(int i = 0; i < size; ++i)
  v[i] = t.v[i];
return *this;
}


Tensor_d1s1::Tensor_d1s1(const Tensor_d1s1 &t) {
size = t.size;
v = new Tensor_d0s1[size];
for(int i = 0; i < size; ++i)
  v[i] = t.v[i];
}













Tensor_d1s2_full::Tensor_d1s2_full(int _size){
size=_size;
v = new Tensor_d0s2[size];
}


Tensor_d1s2_full::Tensor_d1s2_full(const Tensor_d1s2_full &t) {
size = t.size;
v = new Tensor_d0s2[size];
for(int i = 0; i < size; ++i)
  v[i] = t.v[i];
}


Tensor_d1s2_full
Tensor_d1s2_full::operator + (const Tensor_d1s2_full &tens)const {
Tensor_d1s2_full t(size);
for (int k = 0; k < size; k++)
  for (int i = 0; i < 9; i++)
    t[k][i]=v[k][i]+tens[k][i];
return t;
}


Tensor_d1s2_full
Tensor_d1s2_full::operator * (double scal){
Tensor_d1s2_full t(size);
for (int k = 0; k < size; k++)
  for (int i = 0; i < 9; i++)
    t[k][i]=v[k][i]*scal;
return t;
}


Tensor_d1s2_full
operator *(double scal, const Tensor_d1s2_full &tens){
int size = tens.size;
Tensor_d1s2_full t(size);
for (int k = 0; k < size; k++)
  for (int i = 0; i < 9; i++)
    t[k][i]=tens[k][i]*scal;
return t;
};

/*
Tensor_d2s2
Tensor_d1s2_full::operator | (const Tensor_d1s2_full &tens) const {
Tensor_d2s2 t(size);
for (int m = 0; m < size; m++)
  for (int n = 0; n < size; n++)
    for (int i = 0; i < 3; i++)
      for (int j = 0; j < 3; j++)
	for (int k = 0; k < 3; k++)
          t[size*m+n][3*i+j]+=v[m][3*i+k]*tens[n][3*k+j];
return t;
}


Tensor_d1s2_full
Tensor_d1s2_full::operator | (const Tensor_d0s2 &tens) const {
Tensor_d1s2_full t(size);
for (int n = 0; n < size; n++)
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++)
      for (int k = 0; k < 3; k++)
        t[n][3*i+j]+=v[n][3*i+k]*tens[3*k+j];
return t;
}
*/


void
Tensor_d1s2_full::splContractInto(const Tensor &b, Tensor *result) const
{
 
const Tensor_d1s2_full &tens = dynamic_cast<const Tensor_d1s2_full &>(b);
  if(&tens)
  {
    Tensor_d2s2 &t = static_cast<Tensor_d2s2 &>(*result);   
    for (int m = 0; m < size; m++)
      for (int n = 0; n < size; n++)
        for (int i = 0; i < 3; i++)
          for (int j = 0; j < 3; j++){
            t[size*m+n][3*i+j]=0.0;
	    for (int k = 0; k < 3; k++)
              t[size*m+n][3*i+j]+=v[m][3*i+k]*tens[n][3*k+j];
                                      }
  }
  else
  {
    const Tensor_d0s2 &tens1 = static_cast<const Tensor_d0s2 &>(b);
    Tensor_d1s2_full &t = static_cast<Tensor_d1s2_full &>(*result);
    for (int n = 0; n < size; n++)
      for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++){
          t[n][3*i+j]=0.0;
          for (int k = 0; k < 3; k++)
            t[n][3*i+j]+=v[n][3*i+k]*tens1[3*k+j];
                                    }
  }
}




Tensor_d0s2 
Tensor_d1s2_full::operator % (const Tensor_d1s0 &tens) const {
Tensor_d0s2 t;
for (int i = 0; i < 3; i++)
  for (int j = 0; j < 3; j++)
    for (int m = 0; m < size; m++)
       t[3*i+j]+=v[m][3*i+j]*tens[m];
return t;
}


Tensor_d1s2_full &
Tensor_d1s2_full::operator=(const Tensor_d1s2_full &t) {
  if(size != t.size) {
    if(v) delete[] v;
    size = t.size;
    v = new Tensor_d0s2[size];
  }
 for(int i = 0; i < size; ++i)
    v[i] = t.v[i];
 return *this;
}


void
Tensor_d1s2_full::getSpaceTranspose(Tensor_d1s2_full &t){

for (int k = 0; k < size; k++)
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++)
      t[k][3*i+j]=v[k][3*j+i];
}


void 
Tensor_d1s2_full::convertToSym(Tensor_d1s2_Ss23 &t){

for (int k = 0; k < size; k++)   
  for (int i = 0; i < 3; i++)
    for (int j = i; j < 3; j++)
      t[k][i*(5-i)/2+j]=v[k][3*i+j];
}     


Tensor_d1s2_Ss23
Tensor_d1s2_full::symPart(){
Tensor_d1s2_Ss23 t(size);
for (int m=0; m<size; m++)
  for (int i = 0; i < 3; i++)
    for (int j = i; j < 3; j++)
       t[m][i*(5-i)/2+j]=(1./2)*(v[m][3*i+j]+v[m][3*j+i]);
return t;
};














Tensor_d1s2_sparse::Tensor_d1s2_sparse(int _size){
size=_size;
v = new Tensor_d0s2[size/3];
}


Tensor_d1s2_sparse::Tensor_d1s2_sparse(const Tensor_d1s2_sparse &t) {
size = t.size;
v = new Tensor_d0s2[size/3];
for(int i = 0; i < (size/3); ++i)
  v[i] = t.v[i];
}


Tensor_d1s2_sparse
Tensor_d1s2_sparse::operator + (const Tensor_d1s2_sparse &tens)const {
Tensor_d1s2_sparse t(size);
for (int k = 0; k < size/3; k++)
  for (int i = 0; i < 9; i++)
    t[k][i]=v[k][i]+tens[k][i];
return t;
}


Tensor_d1s2_sparse &
Tensor_d1s2_sparse::operator=(const Tensor_d1s2_sparse &t) {
if(size != t.size) {
  if(v) delete[] v;
  size = t.size;
  v = new Tensor_d0s2[size/3];
                   }
//fprintf(stderr, "Un tiens vaut mieux que deux tu l`auras ter\n"); 
for(int i = 0; i < size/3; ++i)
  v[i] = t.v[i];
return *this;
}


Tensor_d0s2 
Tensor_d1s2_sparse::operator % (const Tensor_d1s0 &tens) const {
Tensor_d0s2 t;

for (int j = 0; j < 3; j++)
  for (int m = 0; m < size; m++)
    t[3*(m%3)+j]+=v[(m-(m%3))/3][3*(m%3)+j]*tens[m];
return t;
}

/*
Tensor_d1s2_sparse
Tensor_d1s2_sparse::operator | (const Tensor_d0s2 &tens) const {
Tensor_d1s2_sparse t(size);
for (int n = 0; n < size; n++)
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++)
       t[(n-(n%3))/3][3*(n%3)+i]+=v[(n-(n%3))/3][3*(n%3)+j]*tens[3*j+i];        
return t;
}
*/

void
Tensor_d1s2_sparse::splContractInto(const Tensor &b, Tensor *result) const
{
    const Tensor_d0s2 &tens = static_cast<const Tensor_d0s2 &>(b);
    Tensor_d1s2_sparse &t = static_cast<Tensor_d1s2_sparse &>(*result);
    //fprintf(stderr, "Un tiens vaut mieux que deux tu l`auras\n");
    for (int n = 0; n < size; n++)
      for (int i = 0; i < 3; i++){
        t[(n-(n%3))/3][3*(n%3)+i]=0.0;
        for (int j = 0; j < 3; j++)
           t[(n-(n%3))/3][3*(n%3)+i]+=v[(n-(n%3))/3][3*(n%3)+j]*tens[3*j+i];                                       }
}


Tensor_d1s2_Ss23 
Tensor_d1s2_sparse::symPart(){
Tensor_d1s2_Ss23 t(size);
for (int n = 0; n < size; n++){
  for (int i = (n%3); i < 3; i++)
    t[n][(n%3)*(5-(n%3))/2+i]+=(1./2)*v[(n-(n%3))/3][3*(n%3)+i];
  for (int j = 0; j < (n%3)+1; j++)
     t[n][j*(5-j)/2+(n%3)]+=(1./2)*v[(n-(n%3))/3][3*(n%3)+j];			      				}        
return t;
}


void 
Tensor_d1s2_sparse::getSymSquare(Tensor_d2s2_Sd12s34_full &t){
for (int m = 0; m < size; m++)
  for (int i = 0; i < 3; i++)
    for (int j = i; j < 3; j++)
      t[m*(2*size-m-1)/2+m][i*(5-i)/2+j] = v[(m-(m%3))/3][3*(m%3)+j]*v[(m-(m%3))/3][3*(m%3)+i];
}


void 
Tensor_d1s2_sparse::getSymSquare(Tensor_d2s2_Sd12s34_sparse &t){
for (int m = 0; m < size; m++)
  for (int n = m; n < size; n++)  
    if(((n-m)%3)==0){
      for (int i = 0; i < 3; i++)
        for (int j = i; j < 3; j++)
          t[(size*(size+3)-(size-m+m%3)*(size-m-m%3+3)+2*(n-m))/6][i*(5-i)/2+j] = v[(m-(m%3))/3][3*(m%3)+j]*v[(n-(n%3))/3][3*(n%3)+i];
                  }	           
}






Tensor_d2s0::Tensor_d2s0(int _size, bool isNull){
size=_size; 
if(isNull) 
  v = 0;
else {
  v = new double[size*size]; 
  for (int i = 0; i < (size*size); i++)
    v[i]=0.;
      }
}


Tensor_d2s0::Tensor_d2s0(const Tensor_d2s0 &t) {
size = t.size;
v = new double[size*size];
for(int i = 0; i < (size*size); i++)
  v[i] = t.v[i];
}



Tensor_d2s0
Tensor_d2s0::operator + (const Tensor_d2s0 &tens) const {
if(tens.v == 0)
  return *this;
if(v == 0)
  return tens;
Tensor_d2s0 t(size);
for (int i = 0; i < size*size; i++)
  t[i]=v[i]+tens[i];
return t;
}


Tensor_d2s0 
operator *(double scal, const Tensor_d2s0 &tens){
int size = tens.getSize();
Tensor_d2s0 t(size);
for(int i = 0; i < size*size; ++i)
  t[i] = tens[i]*scal;
return t;
};


Tensor_d2s0 &
Tensor_d2s0::operator=(const Tensor_d2s0 &t) {
if(size != t.size) {
  if(v) delete[] v;
  size = t.size;
  v = new double[size*size];
                    }
for(int i = 0; i < size*size; ++i)
  v[i] = t.v[i];
return *this;
}













/*
Tensor_d2s0_null::Tensor_d2s0_null(const Tensor_d2s0_null &t) {
  size = t.getSize();
  v = 0;
}


*/









Tensor_d2s2::Tensor_d2s2(int _size){size=_size; v = new Tensor_d0s2[size*size];
}


Tensor_d2s2::Tensor_d2s2(const Tensor_d2s2 &t) {
size = t.size;
v = new Tensor_d0s2[size*size];
for(int i = 0; i < size*size; ++i)
  v[i] = t.v[i];
}


Tensor_d2s2
Tensor_d2s2::operator + (const Tensor_d2s2 &tens) const {
Tensor_d2s2 t(size);
for (int k = 0; k < size*size; k++)
  for (int i = 0; i < 9; i++)
    t[k][i]=v[k][i]+tens[k][i];
return t;
}


Tensor_d2s2
Tensor_d2s2::operator * (double scal){
Tensor_d2s2 t(size);
for (int k = 0; k < size*size; k++)
  for (int i = 0; i < 9; i++)
    t[k][i]=v[k][i]*scal;
return t;
}


Tensor_d2s2 
operator *(double scal, const Tensor_d2s2 &tens){
int size=tens.size;
Tensor_d2s2 t(size);
for (int k = 0; k < size*size; k++)
  for (int i = 0; i < 9; i++)
      t[k][i]=tens[k][i]*scal;
return t;
};


Tensor_d2s2 &
Tensor_d2s2::operator=(const Tensor_d2s2 &t) {
if(size != t.size) {
  if(v) delete[] v;
  size = t.size;
  v = new Tensor_d0s2[size*size];
                    }
for(int i = 0; i < size*size; ++i)
  v[i] = t.v[i];
return *this;
}


void
Tensor_d2s2::getSpaceTranspose(Tensor_d2s2 &t){
for (int k = 0; k < size*size; k++)
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++)
      t[k][3*i+j]=v[k][3*j+i];
}


void 
Tensor_d2s2::convertToSym(Tensor_d2s2_Sd12s34_full &t){
for (int i = 0; i < size; i++)
  for (int j = i; j < size; j++)
    for (int k = 0; k < 3; k++)
      for (int l = k; l < 3; l++)
        t[i*(2*size-i-1)/2+j][k*(5-k)/2+l]=v[i*size+j][3*k+l];
}     
















Tensor_d2s2_Sd12s34_null::Tensor_d2s2_Sd12s34_null(const Tensor_d2s2_Sd12s34_null &t) {
size = t.size;
v = new Tensor_d0s2_Ss12[size*(size+1)/2];
}


Tensor_d2s2_Sd12s34_null &
Tensor_d2s2_Sd12s34_null::operator=(const Tensor_d2s2_Sd12s34_null &t) {
if(size != t.size) {
  if(v) delete[] v;
  size = t.size;
  v = new Tensor_d0s2_Ss12[size*(size+1)/2];
                    }
for(int i = 0; i < size*(size+1)/2; ++i)
  v[i] = t.v[i];
return *this;
}













Tensor_d2s2_Sd12s34_full::Tensor_d2s2_Sd12s34_full(int _size){
size=_size;
v = new Tensor_d0s2_Ss12[size*(size+1)/2];
}


Tensor_d2s2_Sd12s34_full::Tensor_d2s2_Sd12s34_full(const Tensor_d2s2_Sd12s34_full &t) {
size = t.size;
v = new Tensor_d0s2_Ss12[size*(size+1)/2];
for(int i = 0; i < size*(size+1)/2; ++i)
  v[i] = t.v[i];
}


Tensor_d2s2_Sd12s34_full &
Tensor_d2s2_Sd12s34_full::operator=(const Tensor_d2s2_Sd12s34_full &t) {
if(size != t.size) {
  if(v) delete[] v;
  size = t.size;
  v = new Tensor_d0s2_Ss12[size*(size+1)/2];
                    }
for(int i = 0; i < size*(size+1)/2; ++i)
  v[i] = t.v[i];
return *this;
}












Tensor_d2s2_Sd12s34_sparse::Tensor_d2s2_Sd12s34_sparse(int _size){
size=_size;
v = new Tensor_d0s2_Ss12[size*(size+3)/6];
}


Tensor_d2s2_Sd12s34_sparse::Tensor_d2s2_Sd12s34_sparse(const Tensor_d2s2_Sd12s34_sparse &t) {
size = t.size;
v = new Tensor_d0s2_Ss12[size*(size+3)/6];
for(int i = 0; i < size*(size+3)/6; ++i)
  v[i] = t.v[i];
}


Tensor_d2s2_Sd12s34_sparse &
Tensor_d2s2_Sd12s34_sparse::operator=(const Tensor_d2s2_Sd12s34_sparse &t) {
if(size != t.size) {
  if(v) delete[] v;
  size = t.size;
  v = new Tensor_d0s2_Ss12[size*(size+3)/6];
                   }
for(int i = 0; i < size*(size+3)/6; ++i)
  v[i] = t.v[i];
return *this;
}

/*
Tensor_d2s0
Tensor_d2s2_Sd12s34_sparse::operator ||(const Tensor_d0s2_Ss12 &tens) const {
Tensor_d2s0 t(size);
for (int i = 0; i < size; i++){
  for (int k = 0; k < 3; k++){
    t[(size+1)*i]+=-tens[k*(5-k)/2+k]*v[i][k*(5-k)/2+k];
    for (int n = k; n < 3; n++)                       
      t[(size+1)*i]+=2*tens[k*(5-k)/2+n]*v[i][k*(5-k)/2+n];
                             };
                              };
return t;
}
*/

void
Tensor_d2s2_Sd12s34_sparse::dblContractInto(const Tensor &b, Tensor *result) const
{
 const Tensor_d0s2_Ss12 &tens = static_cast<const Tensor_d0s2_Ss12 &>(b);
 Tensor_d2s0 &t = static_cast<Tensor_d2s0 &>(*result);
 int i,j;
 
 for (i = 0; i < size; i++)
   for (j = i; j < size; j++){
     t[i*size+j]=0;
     if(((j-i)%3)==0){
       for (int k = 0; k < 3; k++){       
         t[i*size+j]+=-tens[k*(5-k)/2+k]*v[(size*(size+3)-(size-i+i%3)*(size-i-i%3+3)+2*(j-i))/6][k*(5-k)/2+k];
         for (int n = k; n < 3; n++)                       
           t[i*size+j]+=2*tens[k*(5-k)/2+n]*v[(size*(size+3)-(size-i+i%3)*(size-i-i%3+3)+2*(j-i))/6][k*(5-k)/2+n];
                                 }
                      }
       t[j*size+i]= t[i*size+j];           
                              }
}









/*
Tensor_d2s0
Tensor_d2s2_Sd12s34_null::operator ||(const Tensor_d0s2_Ss12 &tens) const {
Tensor_d2s0 t(size,true);
return t;
}
*/
void
Tensor_d2s2_Sd12s34_null::dblContractInto(const Tensor &b, Tensor *result) const
{
 Tensor_d2s0 &t = static_cast<Tensor_d2s0 &>(*result);
 for(int i = 0; i < size; ++i)
   for(int j = 0; j < size; ++j)
     t[i*size+j] = 0.0;
}








/*
Tensor_d2s0
Tensor_d2s2_Sd12s34_full::operator ||(const Tensor_d0s2_Ss12 &tens) const {
Tensor_d2s0 t(size);
for (int i = 0; i < size; i++)
  for (int j = i; j < size; j++){
    for (int k = 0; k < 3; k++){
      t[j+size*i]+=-tens[k*(5-k)/2+k]*v[i*(2*size-i-1)/2+j][k*(5-k)/2+k];
      for (int n = k; n < 3; n++)                       
        t[j+size*i]+=2*tens[k*(5-k)/2+n]*v[i*(2*size-i-1)/2+j][k*(5-k)/2+n];
                               };
                                 };
for (int l = 1; l < size; l++)
  for (int m = 0; m < l; m++)
    t[size*l+m]=t[size*m+l];
return t;
}
*/
void
Tensor_d2s2_Sd12s34_full::dblContractInto(const Tensor &b, Tensor *result) const
{
 const Tensor_d0s2_Ss12 &tens = static_cast<const Tensor_d0s2_Ss12 &>(b);
 Tensor_d2s0 &t = static_cast<Tensor_d2s0 &>(*result);
 for (int i = 0; i < size; i++)
  for (int j = i; j < size; j++){
    t[j+size*i]=0.0;
    for (int k = 0; k < 3; k++){
      t[j+size*i]+=-tens[k*(5-k)/2+k]*v[i*(2*size-i-1)/2+j][k*(5-k)/2+k];
      for (int n = k; n < 3; n++)                       
        t[j+size*i]+=2*tens[k*(5-k)/2+n]*v[i*(2*size-i-1)/2+j][k*(5-k)/2+n];
                               };
                                 };
for (int l = 1; l < size; l++)
  for (int m = 0; m < l; m++)
    t[size*l+m]=t[size*m+l];
}








/*
Tensor_d1s0 
Tensor_d0s2_Ss12::operator ||(const Tensor_d1s2_Ss23 &tens) const{
int size = tens.getSize();
Tensor_d1s0 t(size);
for (int i = 0; i < size; i++)
   t[i]=v[0]*tens[i][0]
       +v[3]*tens[i][3]
       +v[5]*tens[i][5]
       +2*v[1]*tens[i][1]
       +2*v[2]*tens[i][2]
       +2*v[4]*tens[i][4];
return t;
}
*/
void
Tensor_d0s2_Ss12::dblContractInto(const Tensor &b, Tensor *result) const
{
 const Tensor_d1s2_Ss23 &tens = static_cast<const Tensor_d1s2_Ss23 &>(b);
 Tensor_d1s0 &t = static_cast<Tensor_d1s0 &>(*result);
 int size = tens.getSize();
 for (int i = 0; i < size; i++)
   t[i]=v[0]*tens[i][0]
       +v[3]*tens[i][3]
       +v[5]*tens[i][5]
       +2*v[1]*tens[i][1]
       +2*v[2]*tens[i][2]
       +2*v[4]*tens[i][4];
}

void
Tensor_d0s2_Ss12::buildTensorOf(double *state){
for(int i = 0; i < 6; ++i)
  v[i]=state[i];
};


Tensor_d0s2_Ss12
Tensor_d0s2_Ss12::operator + (const Tensor_d0s2_Ss12 &tens) const{
Tensor_d0s2_Ss12 t;
for (int i = 0; i < 6; i++)
  t.v[i]=v[i]+tens[i];
return t;
}


Tensor_d0s2_Ss12
Tensor_d0s2_Ss12::operator - (const Tensor_d0s2_Ss12 &tens) const{
Tensor_d0s2_Ss12 t;
for (int i = 0; i < 6; i++)
  t.v[i]=v[i]-tens[i];
return t;
}


Tensor_d0s2_Ss12 &
Tensor_d0s2_Ss12::operator = (const Tensor_d0s2_Ss12 &t) { 
 for(int i = 0; i < 6; ++i)
    v[i] = t.v[i];
 return *this;
}


Tensor_d0s2_Ss12
operator *(double scal, const Tensor_d0s2_Ss12 &tens){
Tensor_d0s2_Ss12 t;
for (int i = 0; i < 6; i++)
  t[i]=tens[i]*scal;
return t;
};


void 
Tensor_d0s2_Ss12::getDeviation(Tensor_d0s2_Ss12 &t){
double tracethird = (1./3)*(v[0]+v[3]+v[5]);
t[0] = v[0]- tracethird;
t[3] = v[3]- tracethird;
t[5] = v[5]- tracethird;
t[1] = v[1];
t[2] = v[2];
t[4] = v[4];
};


double 
Tensor_d0s2_Ss12::secondInvariant(){
return ((v[0]*v[0]+v[3]*v[3]+v[5]*v[5])/2.+(v[1]*v[1]+v[2]*v[2]+v[4]*v[4]));
};


double 
Tensor_d0s2_Ss12::getTrace(){
return v[0]+v[3]+v[5];
}












/*
Tensor_d1s2_Ss23
Tensor_d0s4_Ss12s34::operator || (const Tensor_d1s2_Ss23 &tens) const {
int size = tens.getSize();
Tensor_d1s2_Ss23 t(size);
for (int m = 0; m < size; m++)
  for (int i = 0; i < 3; i++)
    for (int j = i; j < 3; j++)      
      for (int k = 0; k < 3; k++){
        t[m][i*(5-i)/2+j]+=-v[(i*(5-i)/2+j)][k*(5-k)/2+k]*tens[m][k*(5-k)/2+k];
        for (int l = k; l < 3; l++)
          t[m][i*(5-i)/2+j]+=2*v[(i*(5-i)/2+j)][k*(5-k)/2+l]*tens[m][k*(5-k)/2+l];
                                 };
return t;
}


Tensor_d0s2_Ss12
Tensor_d0s4_Ss12s34::operator || (const Tensor_d0s2_Ss12 &tens) const {
Tensor_d0s2_Ss12 t;
for (int i = 0; i < 3; i++)
  for (int j = i; j < 3; j++)      
    for (int k = 0; k < 3; k++){
      for (int l = k; l < 3; l++)
        t[i*(5-i)/2+j]+=v[(i*(5-i)/2+j)][k*(5-k)/2+l]*tens[k*(5-k)/2+l];
                               };      
return t;
}
*/

void
Tensor_d0s4_Ss12s34::dblContractInto(const Tensor &b, Tensor *result) const
{
  const Tensor_d0s2_Ss12 *tens = dynamic_cast<const Tensor_d0s2_Ss12 *>( &b);
  if(tens)
  {
    Tensor_d0s2_Ss12 &t = static_cast<Tensor_d0s2_Ss12 &>(*result);   
    for (int i = 0; i < 3; i++)
      for (int j = i; j < 3; j++){
        t[i*(5-i)/2+j]=0.0;     
        for (int k = 0; k < 3; k++){
          for (int l = k; l < 3; l++)
            t[i*(5-i)/2+j]+=v[(i*(5-i)/2+j)][k*(5-k)/2+l]*(*tens)[k*(5-k)/2+l];
                                   } 
                                  }     
  }
  else
  {
    const Tensor_d1s2_Ss23 &tens = static_cast<const Tensor_d1s2_Ss23 &>(b);
    Tensor_d1s2_Ss23 &t = static_cast<Tensor_d1s2_Ss23 &>(*result);
    int size = tens.getSize();
    for (int m = 0; m < size; m++)
      for (int i = 0; i < 3; i++)
        for (int j = i; j < 3; j++){
          t[m][i*(5-i)/2+j]=0.0;     
          for (int k = 0; k < 3; k++)                           
            for (int l = k; l < 3; l++)              t[m][i*(5-i)/2+j]+=v[(i*(5-i)/2+j)][k*(5-k)/2+l]*tens[m][k*(5-k)/2+l];
                                    }
  }



}



void
Tensor_d0s4_Ss12s34::reOrder(Tensor_d0s4_Ss12s34 tm){

for (int i = 0; i<6; ++i){

  v[i][0]=tm[i][0];
  v[i][1]=2*tm[i][1];
  v[i][2]=2*tm[i][2];
  v[i][3]=tm[i][3];
  v[i][4]=2*tm[i][4];
  v[i][5]=tm[i][5];
                        }


}




















Tensor_d1s2_Ss23 
Tensor_d1s2_Ss23::operator + (const Tensor_d1s2_Ss23 &tens) const{
Tensor_d1s2_Ss23 t(size);
for (int k = 0; k < size; k++)
  for (int i = 0; i < 6; i++)
      t.v[k][i]=v[k][i]+tens[k][i];
return t;
}


Tensor_d1s2_Ss23 &
Tensor_d1s2_Ss23::operator=(const Tensor_d1s2_Ss23 &t) {
if(size != t.size) {
  if(v) delete[] v;
  size = t.size;
  v = new Tensor_d0s2_Ss12[size];
                    }
//fprintf(stderr, "Un tiens vaut mieux que deux tu l`auras ter\n");
for(int i = 0; i < size; ++i)
  v[i] = t.v[i];
return *this;
}


Tensor_d1s2_Ss23::Tensor_d1s2_Ss23(int _size){
size=_size;
v = new Tensor_d0s2_Ss12[size];
}


Tensor_d1s2_Ss23::Tensor_d1s2_Ss23(const Tensor_d1s2_Ss23 &t) {
size = t.getSize();
v = new Tensor_d0s2_Ss12[size];
for(int i = 0; i < size; ++i)
  v[i] = t.v[i];
}

/*
Tensor_d2s0
Tensor_d1s2_Ss23::operator ||(const Tensor_d1s2_Ss23 &tens) const{
Tensor_d2s0 t(size);
for (int m = 0; m < size; m++)
  for (int n = m; n < size; n++)
    for (int i = 0; i < 3; i++){
      t[n+size*m]+=-v[m][i*(5-i)/2+i]*tens[n][i*(5-i)/2+i];
      for (int j = i; j < 3; j++)
        t[n+size*m]+=2*v[m][i*(5-i)/2+j]*tens[n][i*(5-i)/2+j];
                                };
for (int o = 1; o < size; o++)
  for (int p = 0; p < o; p++)
    t[p+size*o]=t[o+size*p];
return t;
}
*/

void
Tensor_d1s2_Ss23::dblContractInto(const Tensor &b, Tensor *result) const
{
 const Tensor_d1s2_Ss23 &tens = static_cast<const Tensor_d1s2_Ss23 &>(b);
 Tensor_d2s0 &t = static_cast<Tensor_d2s0 &>(*result);
 for (int m = 0; m < size; m++)
  for (int n = m; n < size; n++){
    t[n+size*m]=0.0;
    for (int i = 0; i < 3; i++){
      t[n+size*m]+=-v[m][i*(5-i)/2+i]*tens[n][i*(5-i)/2+i];
      for (int j = i; j < 3; j++)
        t[n+size*m]+=2*v[m][i*(5-i)/2+j]*tens[n][i*(5-i)/2+j];
                               }
                                 } 
for (int o = 1; o < size; o++)
  for (int p = 0; p < o; p++)
    t[p+size*o]=t[o+size*p];
}
