#ifndef _NODE_SPACE_H_
#define _NODE_SPACE_H_
#include <stdio.h>
#include <stdlib.h>

class DoubleContraction;
class Contraction;

class Tensor {
  public:
    Contraction operator | (Tensor &);
    DoubleContraction operator || (Tensor &);
    //virtual Tensor operator + 
    virtual void dblContractInto(const Tensor &, Tensor *) const ;
    virtual void splContractInto(const Tensor &, Tensor *) const ;
    Tensor &operator= (const DoubleContraction &);
    Tensor &operator= (const Contraction &);
};


 
class  Contraction{

   Tensor &a, &b;

 public:

   Contraction(Tensor &_a,Tensor &_b) : a(_a), b(_b) {};
   void assignTo(Tensor *result) const { a.splContractInto(b, result); }
};    


class DoubleContraction{

   Tensor &a, &b;
 public:

   DoubleContraction(Tensor &_a,Tensor &_b) : a(_a), b(_b) {};
   void assignTo(Tensor *result) const { a.dblContractInto(b, result); }

};


inline
Tensor &Tensor::operator= (const DoubleContraction &dc)
{
 dc.assignTo(this);
 return *this;
}


inline
Tensor &Tensor::operator= (const Contraction &dc)
{
 dc.assignTo(this);
 return *this;
}

 


class Tensor_d0s2_Ss12 ;

class Tensor_d0s4_Ss12s34 ;

class Tensor_d2s0_Sd12 ;

class Tensor_d2s0_null;

class Tensor_d2s2_Sd12s34 ;

class Tensor_d2s2_Sd12s34_full ;

class Tensor_d2s2_Sd12s34_sparse ;

class Tensor_d2s2_Sd12s34_null;

class Tensor_d1s1;

class Tensor_d1s2;

class Tensor_d1s2_Ss23 ;

class Tensor_d1s2_full;

class Tensor_d1s2_sparse;

class Tensor_d2s2;






class Tensor_d0s1 : public Tensor {

protected:

double v[3];

public:
    Tensor_d0s1(){for (int i = 0; i<3; i++) v[i]=0;}
    double &operator[] (int i){return v[i];}
    double operator[] (int i) const {return v[i];}
    double &operator() (int i){return v[i];}
    double operator() (int i) const {return v[i];}
};
 



class Tensor_d0s2 : public Tensor {

protected:

double v[9];

public:
    Tensor_d0s2(){for (int i = 0; i<9; i++) v[i]=0;}
    Tensor_d0s2 operator + (const Tensor_d0s2 &) const;
    Tensor_d0s2 operator * (double scal);
    //Tensor_d0s2 operator | (const Tensor_d0s2 &)const;
    //Tensor_d1s2_full operator | (const Tensor_d1s2_full &) const;
    //Tensor_d1s2_full operator | (const Tensor_d1s2_sparse &) const;
    Tensor_d0s2 &operator = (const Contraction &dc);
    void splContractInto(const Tensor &b, Tensor *result) const;
    double &operator[] (int i){return v[i];}
    double operator[] (int i) const {return v[i];}
    double &operator() (int i, int j){return v[3*i+j];}
    double operator() (int i, int j) const{return v[3*i+j];}
    void getDeterminant(double &det);
    void getInverse(Tensor_d0s2 &t);
    void getTranspose(Tensor_d0s2 &t);
    void convertToSym(Tensor_d0s2_Ss12 &t);
    friend Tensor_d0s2 operator *(double d, const Tensor_d0s2 &t);
};
 
Tensor_d0s2 operator *(double, const Tensor_d0s2 &);


inline
Tensor_d0s2 &Tensor_d0s2::operator= (const Contraction &dc)
{
 dc.assignTo(this);
 return *this;
}



class Tensor_d0s4 : public Tensor {
public:
    Tensor_d0s4 operator + (Tensor_d0s4 &);
    Tensor_d0s4 operator + (const Tensor_d0s4 &) const;
};

class Tensor_d1s0 : public Tensor {

protected:

int size;
double *v;

public:
    Tensor_d1s0(){size = 0; v = 0;};
    Tensor_d1s0(int _size);
    Tensor_d1s0(const Tensor_d1s0 &t);
    double &operator[] (int i) {return v[i];}
    double operator[] (int i) const {return v[i];}
    double &operator() (int i) {return v[i];}
    double operator() (int i) const {return v[i];}
    Tensor_d1s0 operator + (const Tensor_d1s0 &) const;
    Tensor_d1s0 &operator=(const Tensor_d1s0 &);
    Tensor_d1s0 &operator=(const DoubleContraction &);
    int getSize() const {return size;}
    friend Tensor_d1s0 operator *(double d, const Tensor_d1s0 &t);
};

Tensor_d1s0 operator *(double d, const Tensor_d1s0 &t);


inline
Tensor_d1s0 &Tensor_d1s0::operator= (const DoubleContraction &dc)
{
 dc.assignTo(this);
 return *this;
}



class Tensor_d1s1 : public Tensor {

protected:

int size;
Tensor_d0s1 *v;

public:
    Tensor_d1s1(){ size = 0; v = 0;};
    Tensor_d1s1(int _size);
    Tensor_d1s1(const Tensor_d1s1 &t);
    ~Tensor_d1s1() { if (v) delete [] v; }
    Tensor_d0s1 &operator[] (int i) {return v[i];}
    Tensor_d0s1 operator[] (int i) const {return v[i];}  
    double operator() (int i, int j) const;
    double &operator() (int i, int j);
    Tensor_d1s1 &operator=(const Tensor_d1s1 &);
    int getSize() const {return size;}
};

inline double
Tensor_d1s1::operator() (int i, int j) const {
#ifdef DEBUG
if(i>=size){throw "first index out of range\n";} 
#endif
return v[i][j];
};

inline double&
Tensor_d1s1::operator() (int i, int j) {
#ifdef DEBUG
if(i>=size){throw "first index out of range\n";} 
#endif
return v[i][j];
};



class Tensor_d1s2 : public Tensor {

protected:

int size;

public:

    virtual Tensor_d1s2_Ss23 symPart()=0;
};





class Tensor_d1s2_full : public Tensor_d1s2 {

protected:

Tensor_d0s2 *v;

public:
    Tensor_d1s2_full(){ size = 0; v = 0;};
    Tensor_d1s2_full(int _size);
    Tensor_d1s2_full(const Tensor_d1s2_full &t);
    ~Tensor_d1s2_full() { if (v) delete [] v; }
    Tensor_d1s2_full operator + (const Tensor_d1s2_full &) const;
    Tensor_d1s2_full operator * (double scal);
    Tensor_d0s2 &operator[] (int i) {return v[i];}
    Tensor_d0s2 operator[] (int i) const {return v[i];}
    double &operator() (int i,int j, int k);
    double operator() (int i,int j, int k) const;
    Tensor_d1s2_full &operator= (const Tensor_d1s2_full &);
    Tensor_d1s2_full &operator= (const Contraction &);
    //Tensor_d2s2 operator | (const Tensor_d1s2_full &) const;
    //Tensor_d1s2_full operator | (const Tensor_d0s2 &) const;
    void splContractInto(const Tensor &b, Tensor *result) const;
    Tensor_d0s2 operator % (const Tensor_d1s0 &) const;
    void getSpaceTranspose(Tensor_d1s2_full &t);
    void convertToSym(Tensor_d1s2_Ss23 &t);
    Tensor_d1s2_Ss23 symPart();
    int getSize() const {return size;}
    friend Tensor_d1s2_full operator *(double d, const Tensor_d1s2_full &t);
};

inline double 
Tensor_d1s2_full::operator()(int i, int j, int k) const {
#ifdef DEBUG
if(i>=size){throw"first index out of range\n";} 
#endif
return v[i][3*j+k];
};

inline double& 
Tensor_d1s2_full::operator()(int i, int j, int k) {
#ifdef DEBUG
if(i>=size){throw"first index out of range\n");} 
#endif
return v[i][3*j+k];
};




Tensor_d1s2_full operator *(double d, const Tensor_d1s2_full &t);


inline
Tensor_d1s2_full &Tensor_d1s2_full::operator= (const Contraction &dc)
{
 dc.assignTo(this);
 return *this;
}




class Tensor_d1s2_sparse : public Tensor_d1s2 {

protected:

Tensor_d0s2 *v;

public:
    Tensor_d1s2_sparse(){ size = 0; v = 0;};
    Tensor_d1s2_sparse(int _size);
    Tensor_d1s2_sparse(const Tensor_d1s2_sparse &t);
    ~Tensor_d1s2_sparse() { if (v) delete [] v; }
    Tensor_d1s2_sparse operator + (const Tensor_d1s2_sparse &) const;
    Tensor_d0s2 &operator[] (int i) {return v[i];}
    Tensor_d0s2 operator[] (int i) const {return v[i];}
    double &operator() (int i, int j, int k);
    double operator() (int i, int j, int k) const;
    Tensor_d1s2_sparse &operator= (const Tensor_d1s2_sparse &);
    Tensor_d1s2_sparse &operator= (const Contraction &);
    //Tensor_d1s2_sparse operator | (const Tensor_d0s2 &) const;
    void splContractInto(const Tensor &b, Tensor *result) const;
    Tensor_d0s2 operator % (const Tensor_d1s0 &) const;
    Tensor_d1s2_Ss23 symPart();
    void getSymSquare(Tensor_d2s2_Sd12s34_full &t);
    void getSymSquare(Tensor_d2s2_Sd12s34_sparse &t);
    int getSize() const {return size;}
    friend Tensor_d1s2_sparse operator *(double d, const Tensor_d1s2_sparse &t);
};

inline double 
Tensor_d1s2_sparse::operator()(int i, int j, int k) const {
#ifdef DEBUG
if(i>=size){throw"first index out of range\n");} 
#endif
if((j-k)%3==0){return v[(i-(i%3))/3][3*(i%3)+k];}
else{return 0;}
};


inline double& 
Tensor_d1s2_sparse::operator()(int i, int j, int k) {
#ifdef DEBUG
if(i>=size){throw"first index out of range\n");} 
#endif
if((j-k)%3==0){return v[(i-(i%3))/3][3*(i%3)+k];}
else{throw"Error : the tensor is sparse ! Check the constructor.\n";}
};





Tensor_d1s2_sparse operator *(double d, const Tensor_d1s2_sparse &t);


inline
Tensor_d1s2_sparse &Tensor_d1s2_sparse::operator= (const Contraction &dc)
{
//fprintf(stderr, "Un tiens vaut mieux que deux tu l`auras\n"); 
 dc.assignTo(this);
 return *this;
}




class Tensor_d2s2 : public Tensor {

protected:

int size;
Tensor_d0s2 *v;

public:
    Tensor_d2s2(){size = 0; v = 0;};
    Tensor_d2s2(int _size);
    Tensor_d2s2(const Tensor_d2s2 &t);
    ~Tensor_d2s2() { if (v) delete [] v; }
    Tensor_d2s2 operator + (const Tensor_d2s2 &) const;
    Tensor_d2s2 operator * (double scal) ;
    Tensor_d0s2 &operator[] (int i) {return v[i];}
    Tensor_d0s2 operator[] (int i) const {return v[i];}
    double &operator() (int i, int j, int k, int l);
    double operator() (int i, int j, int k, int l) const;
    Tensor_d2s2 &operator=(const Tensor_d2s2 &);
    int getSize() const {return size;}
    void getSpaceTranspose(Tensor_d2s2 &t);
    void convertToSym(Tensor_d2s2_Sd12s34_full &t);
    friend Tensor_d2s2 operator *(double d, const Tensor_d2s2 &t);
};



inline double 
Tensor_d2s2::operator()(int i, int j, int k, int l) const {
#ifdef DEBUG
if(i>=size||j>=size){throw"Index out of range\n";} 
#endif
return v[i*size+j][k*3+l];
};


inline double& 
Tensor_d2s2::operator()(int i, int j, int k, int l) {
#ifdef DEBUG
if(i>=size||j>=size){throw"Index out of range\n";} 
#endif
return v[i*size+j][k*3+l];
};



Tensor_d2s2 operator *(double d, const Tensor_d2s2 &t);






class Tensor_d2s0 : public Tensor {

protected:
  int size;
  double *v;

public:
    Tensor_d2s0(){size = 0; v = 0;}
    Tensor_d2s0(int _size, bool isNull );
    Tensor_d2s0(const Tensor_d2s0 &t);
    ~Tensor_d2s0() { if (v) delete [] v; }
    Tensor_d2s0 operator + (const Tensor_d2s0 &) const;
    double &operator[] (int i) {return v[i];}
    double operator[] (int i) const {return v[i];}
    double &operator() (int i, int j);
    double operator() (int i, int j) const;
    Tensor_d2s0 &operator=(const Tensor_d2s0 &);
    Tensor_d2s0 &operator=(const DoubleContraction &);
    int getSize() const {return size;}
    friend Tensor_d2s0 operator *(double d, const Tensor_d2s0 &t);
};


inline double 
Tensor_d2s0::operator()(int i, int j) const {
#ifdef DEBUG
if(i>=size||j>=size){throw"Index out of range\n";} 
#endif
return v[i*size+j];
};

inline double& 
Tensor_d2s0::operator()(int i, int j) {
#ifdef DEBUG
if(i>=size||j>=size){throw"Index out of range\n");} 
#endif
return v[i*size+j];
};


Tensor_d2s0 operator *(double, const Tensor_d2s0 &);

inline
Tensor_d2s0 &Tensor_d2s0::operator= (const DoubleContraction &dc)
{
 dc.assignTo(this);
 return *this;
}


/*
class Tensor_d2s0_null : public Tensor_d2s0 {

protected:
  int size;
  double *v;

public:
    Tensor_d2s0_null(){size = 0; v=0;}
    Tensor_d2s0_null(const Tensor_d2s0_null &t);
    ~Tensor_d2s0_null() { if (v) delete [] v; }
    Tensor_d2s0_null(int _size){size=_size; v = 0;}
    int getSize() const {return size;}
};
*/






class Tensor_d0s2_Ss12 : public Tensor {//Stress ans Strain  

protected:

double v[6]; 

public:
      Tensor_d0s2_Ss12(){for (int i = 0; i<6; i++) v[i]=0;}
      Tensor_d1s0 operator ||(const Tensor_d1s2_Ss23 &) const;
      //Tensor_d2s0 operator ||(const Tensor_d2s2_Sd12s34_full &) const;
      //Tensor_d2s0 operator ||(const Tensor_d2s2_Sd12s34_sparse &) const;
      double &operator[] (int i){return v[i];}
      double operator[] (int i) const {return v[i];}
      double &operator() (int i, int j);
      double operator() (int i, int j) const;
      Tensor_d0s2_Ss12 operator + (const Tensor_d0s2_Ss12  &) const;
      Tensor_d0s2_Ss12 operator - (const Tensor_d0s2_Ss12  &) const;
      Tensor_d0s2_Ss12 &operator = (const Tensor_d0s2_Ss12  &);
      Tensor_d0s2_Ss12 &operator = (const DoubleContraction &);
      void buildTensorOf(double *state);
      void getDeviation(Tensor_d0s2_Ss12 &t);
      double getTrace();
      double secondInvariant();

      void dblContractInto(const Tensor &, Tensor *) const;
      friend Tensor_d0s2_Ss12 operator *(double d, const Tensor_d0s2_Ss12 &t);
};

inline double 
Tensor_d0s2_Ss12::operator()(int i, int j) const {
#ifdef DEBUG
if(j>=i){return v[i*(5-i)/2+j];}
#endif
return v[j*(5-j)/2+i];
};

inline double&
Tensor_d0s2_Ss12::operator()(int i, int j){
#ifdef DEBUG
if(i>j){throw"Error : symetric tensor indices are triangular\n";}
#endif
return v[i*(5-i)/2+j];
};




Tensor_d0s2_Ss12 operator *(double, const Tensor_d0s2_Ss12 &);


inline
Tensor_d0s2_Ss12 &Tensor_d0s2_Ss12::operator= (const DoubleContraction &dc)
{
 dc.assignTo(this);
 return *this;
}



class Tensor_d0s4_Ss12s34 : public Tensor{//Tangent material

protected:  

double v[6][6];

public:
     Tensor_d0s4_Ss12s34(){for (int i = 0; i<6; i++) 
                             for(int j = 0; j<6; j++) v[i][j]=0;}
     //Tensor_d1s2_Ss23  operator || (const Tensor_d1s2_Ss23 &) const;
     //Tensor_d0s2_Ss12  operator || (const Tensor_d0s2_Ss12 &) const;
     void dblContractInto(const Tensor &, Tensor *) const;
     void reOrder(Tensor_d0s4_Ss12s34 tm);
     const double *operator[] (int i) const {return v[i];}
     double *operator[] (int i) {return v[i];}
     double &operator() (int i, int j, int k, int l);
     double operator() (int i, int j, int k, int l) const;
};

inline double 
Tensor_d0s4_Ss12s34::operator()(int i, int j, int k, int l) const {
if(j>=i){
  if(l>=k){
    return v[i*(5-i)/2+j][k*(5-k)/2+l];
          }
  else{
    return v[i*(5-i)/2+j][l*(5-l)/2+k];
      }
         }
else{
   if(l>=k){
    return v[j*(5-j)/2+i][k*(5-k)/2+l];
           }
  else{
    return v[j*(5-j)/2+i][l*(5-l)/2+k];
      }
    }
};

inline double& 
Tensor_d0s4_Ss12s34::operator()(int i, int j, int k, int l){
if((j>=i)&&(l>=k)){
  return v[i*(5-i)/2+j][k*(5-k)/2+l];
                  }
else{throw"Error : symetric tensor indices are triangular\n";}

};









class Tensor_d2s2_Sd12s34 : public Tensor {//dB

protected:  

int size;

public:

    Tensor_d2s2_Sd12s34(){};
    Tensor_d2s2_Sd12s34(int _size){size = _size;}
    virtual int getSize(){return size;}
    //virtual Tensor_d2s0 operator ||(const Tensor_d0s2_Ss12 &) const = 0;
    virtual void dblContractInto(const Tensor &, Tensor *) const = 0;   
};








class Tensor_d2s2_Sd12s34_full : public Tensor_d2s2_Sd12s34 {//dB

protected:  

Tensor_d0s2_Ss12 *v;

public:

    Tensor_d2s2_Sd12s34_full(){size = 0; v=0;};
    Tensor_d2s2_Sd12s34_full(int _size);
    Tensor_d2s2_Sd12s34_full(const Tensor_d2s2_Sd12s34_full &t);
    ~Tensor_d2s2_Sd12s34_full() { if (v) delete [] v; }
    Tensor_d0s2_Ss12 &operator[] (int i) {return v[i];}
    Tensor_d0s2_Ss12 operator[] (int i) const {return v[i];}
    double &operator() (int i, int j, int k, int l);
    double operator() (int i, int j, int k, int l) const;
    Tensor_d2s2_Sd12s34_full &operator=(const Tensor_d2s2_Sd12s34_full&);
    //Tensor_d2s0 operator ||(const Tensor_d0s2_Ss12 &) const;
    void dblContractInto(const Tensor &, Tensor *) const;
    int getSize(){return size;}
    int getSize() const {return size;}

};

inline double 
Tensor_d2s2_Sd12s34_full::operator()(int i, int j, int k, int l) const {
#ifdef DEBUG
if(i>=size||j>=size){throw"Index out of range\n";} 
#endif
if(j>=i){
  if(l>=k){
    return v[i*(2*size-i-1)/2+j][k*(5-k)/2+l];
          }
  else{
    return v[i*(2*size-i-1)/2+j][l*(5-l)/2+k];
      }
        }
else{
   if(l>=k){
    return v[j*(2*size-j-1)/2+i][k*(5-k)/2+l];
           }
  else{
    return v[j*(2*size-j-1)/2+i][l*(5-l)/2+k];
      }
    }
};

inline double& 
Tensor_d2s2_Sd12s34_full::operator()(int i, int j, int k, int l){
#ifdef DEBUG
if(i>=size||j>=size){throw"Index out of range\n";} 
#endif
if((j>=i)&&(l>=k)){
  return v[i*(2*size-i-1)/2+j][k*(5-k)/2+l];
                  }
else{throw"Error : symetric tensor indices are triangular\n";}

};



class Tensor_d2s2_Sd12s34_sparse : public Tensor_d2s2_Sd12s34 {//dB

protected:  

Tensor_d0s2_Ss12 *v;

public:

    Tensor_d2s2_Sd12s34_sparse(){size = 0; v=0;};
    Tensor_d2s2_Sd12s34_sparse(int _size);
    Tensor_d2s2_Sd12s34_sparse(const Tensor_d2s2_Sd12s34_sparse &t);
    ~Tensor_d2s2_Sd12s34_sparse() { if (v) delete [] v; }
    Tensor_d0s2_Ss12 &operator[] (int i) {return v[i];}
    Tensor_d0s2_Ss12 operator[] (int i) const {return v[i];}
    double &operator() (int i, int j, int k, int l);
    double operator() (int i, int j, int k, int l) const;
    Tensor_d2s2_Sd12s34_sparse &operator=(const Tensor_d2s2_Sd12s34_sparse &);
    //Tensor_d2s0 operator ||(const Tensor_d0s2_Ss12 &) const;
    void dblContractInto(const Tensor &, Tensor *) const;
    int getSize(){return size;}
    int getSize() const {return size;}

};

inline double 
Tensor_d2s2_Sd12s34_sparse::operator()(int m, int n, int i, int j) const {
#ifdef DEBUG
if(m>=size||n>=size){throw"Index out of range\n";} 
#endif
if((m-n)%3==0){
  if(n>=m){
    if(j>=i){    
      return v[(size*(size+3)-(size-m+m%3)*(size-m-m%3+3)+2*(n-m))/6][i*(5-i)/2+j];
            }                   
    else{
      return v[(size*(size+3)-(size-m+m%3)*(size-m-m%3+3)+2*(n-m))/6][j*(5-j)/2+i];
        }
          }
  else{
     if(j>=i){
      return v[(size*(size+3)-(size-n+n%3)*(size-n-n%3+3)+2*(m-n))/6][i*(5-i)/2+j];
             }
    else{
      return v[(size*(size+3)-(size-n+n%3)*(size-n-n%3+3)+2*(m-n))/6][j*(5-j)/2+i];
        }
      }
              }
else{return 0;}
};

inline double& 
Tensor_d2s2_Sd12s34_sparse::operator()(int m, int n, int i, int j){
#ifdef DEBUG
if(m>=size||n>=size){throw"Index out of range\n";} 
#endif
if((m-n)%3==0){
  if((n>=m)&&(j>=i)){   
      return v[(size*(size+3)-(size-m+m%3)*(size-m-m%3+3)+2*(n-m))/6][i*(5-i)/2+j];
                    }                   
  else{throw"Error : symetric tensor indices are triangular\n";}
              }
else{throw"Error : the tensor is sparse ! Check the constructor.\n";}
            
};



class Tensor_d2s2_Sd12s34_null : public Tensor_d2s2_Sd12s34 {

protected: 

Tensor_d0s2_Ss12 *v;

public:
    Tensor_d2s2_Sd12s34_null(){size = 0; v=0;};
    Tensor_d2s2_Sd12s34_null(int _size){size=_size; v=0;}
    Tensor_d2s2_Sd12s34_null(const Tensor_d2s2_Sd12s34_null &t); 
    ~Tensor_d2s2_Sd12s34_null() { if (v) delete [] v; }
    Tensor_d2s2_Sd12s34_null &operator=(const Tensor_d2s2_Sd12s34_null &);
    //double operator()(int i, int j, int k, int l) const;
    //double &operator()(int i, int j, int k, int l) const;
    //Tensor_d2s0 operator ||(const Tensor_d0s2_Ss12 &) const;
    void dblContractInto(const Tensor &, Tensor *) const;
    int getSize(){return size;}
    int getSize() const {return size;}
    
};


/*
inline double 
Tensor_d2s2_Sd12s34_null::operator()(int m, int n, int i, int j)const {
#ifdef DEBUG
if(m>=size||n>=size){fprintf(stderr, "Index out of range\n");exit(10);} 
#endif
return 0;
};

inline double& 
Tensor_d2s2_Sd12s34_null::operator()(int m, int n, int i, int j) const {
#ifdef DEBUG
if(m>=size||n>=size){fprintf(stderr, "Index out of range\n");exit(10);} 
#endif
fprintf(stderr, "Error : A null tensor doesn`t need to be filled !.\n");
exit(10);   
};
*/


class Tensor_d1s2_Ss23 : public Tensor {//B

protected:

int size;
Tensor_d0s2_Ss12 *v;

   public:
      Tensor_d1s2_Ss23(){};
      Tensor_d1s2_Ss23(int _size);
      Tensor_d1s2_Ss23(const Tensor_d1s2_Ss23 &t);
      ~Tensor_d1s2_Ss23() { if (v) delete [] v; }
      //Tensor_d2s0 operator ||(const Tensor_d1s2_Ss23 &) const;
      void dblContractInto(const Tensor &, Tensor *) const;
      Tensor_d0s2_Ss12 &operator[] (int i) {return v[i];}
      Tensor_d0s2_Ss12 operator[] (int i) const {return v[i];}
      double &operator() (int i, int j, int k);
      double operator() (int i, int j, int k) const;
      Tensor_d1s2_Ss23 operator + (const Tensor_d1s2_Ss23 &) const;
      Tensor_d1s2_Ss23 &operator=(const Tensor_d1s2_Ss23&);
      int getSize() const {return size;}
};


inline double 
Tensor_d1s2_Ss23::operator()(int i, int j, int k) const {
#ifdef DEBUG
if(i>=size){throw"Index out of range\n");} 
#endif
if(k>=j){return v[i][j*(5-j)/2+k];}
else{return v[i][k*(5-k)/2+j];}
};

inline double&
Tensor_d1s2_Ss23::operator()(int i, int j, int k){
#ifdef DEBUG
if(i>=size){throw"Index out of range\n");} 
#endif
if(j>k){throw"Error : symetric tensor indices are triangular\n";}
return v[i][j*(5-j)/2+k];
};




#endif
