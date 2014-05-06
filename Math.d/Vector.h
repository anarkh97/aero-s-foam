#ifndef _VECTOR_H_
#define _VECTOR_H_
#include <Utils.d/NodeSpaceArray.h>
#include <Utils.d/MyComplex.h>
#include <iostream>

class GeomState;
class CoordSet;
class SingleInfo;

struct BlockInfo {
  int blocklen;
  int numblocks;
};


class SingleInfo {
 
 //CD: for Pita
 public:
  int len;

  SingleInfo(int l = 0)  { len = l; }

  int totLen() const { return len; }

};

// StackVector should be the base class without a destructor and
// Vector derived from it
template<class Scalar>
class GenVector {
  protected:
   int   len;
   Scalar *d;
   bool myMemory;
   int n; // block size
  public:
   // Constructors
   GenVector() { len = 0; d = 0; myMemory = true; }
   GenVector(int length); 
   GenVector(const SingleInfo &inf);       //CD: for Pita
   GenVector(const GenVector<Scalar> &);
   GenVector(const GenVector<Scalar> &, Scalar initialValue);
   GenVector(int length, Scalar initialValue);
   GenVector(Scalar  *v, int length, bool myMemory = true);
   GenVector(int length, Scalar *v, bool myMemory = true); //CBM
   GenVector(int size, GeomState *, CoordSet*);
   GenVector(const GenVector<Scalar> &, int nr, int sr);
   GenVector(const GenVector<Scalar> &, int nr, int *rows);
   GenVector(const BlockInfo &inf);

   virtual ~GenVector(); //HB: made the destructor virtual to deal with
                         //    the derived class GenStackVector

   void setmyMemory(bool _myMemory) {myMemory = _myMemory;}

   int size() const { return len; }
   int info() const { return len; } // JC: added to have same behaviour as NewVec::Vec in FVector.h 

   void initialize(int length);

   void copy(const GenVector<Scalar> &v);
   void copy(const Scalar *v);
   void copy(const Scalar  value);

   Scalar operator * (const GenVector<Scalar> &) const;
   Scalar operator ^ (const GenVector<Scalar> &);
   GenVector operator / (const GenVector<Scalar> &);
   GenVector operator + (const GenVector<Scalar> &);
   GenVector operator - (const GenVector<Scalar> &);
   Scalar &operator[] (int i) const;

   GenVector &operator=(const GenVector<Scalar> &);     // v2 = v1;
   GenVector &operator=(const Scalar c);     // v1 = 0.0;
   GenVector &operator=(const Scalar *data); // v1 = data;

   void operator*=(const Scalar c);
   void operator*(const Scalar c);
   void operator+=(const GenVector<Scalar> &v2);
   void operator-=(const GenVector<Scalar> &v2);
   void operator/=(const GenVector<Scalar> &v2);

   void linAdd(Scalar alpha, const GenVector<Scalar> &v);
   void linAdd(Scalar alpha, const GenVector<Scalar> &, Scalar beta,
               const GenVector<Scalar> &);
   void linC(Scalar c, const GenVector<Scalar> &v);
   void linC(const GenVector<Scalar> &v, Scalar c);
   void linC(const GenVector<Scalar> &, Scalar, const GenVector<Scalar> & );
   void linC(Scalar, const GenVector<Scalar> &, Scalar, const GenVector<Scalar> &);
   void linC(Scalar, const GenVector<Scalar> &, Scalar, const GenVector<Scalar> &, Scalar, const GenVector<Scalar> &);
   void swap(GenVector<Scalar> &);
   Scalar* data() const { return d; }

   void updateBlock(int ii, Scalar c, GenVector<Scalar> &v);
   void copyBlock(GenVector<Scalar> &v, int ii);
   void addBlockSqr(int ii, Scalar c, GenVector<Scalar> &v);
   void computeSqrt();
   void computeRealz(int ii, Scalar c, GenVector<Scalar> &v);

   void print(const char *msg = "",const char *msg2="v");
   GenVector cross(GenVector<Scalar> &v2);
   double magnitude();
   double absMax();
   void zero();
   void clean_up();
   void mult(GenVector<Scalar> &v);
   void diff(GenVector<Scalar> &v1, GenVector<Scalar> &v2);
   void add(GenVector<Scalar>&, int);
   void add(Scalar *v);
   void add(GenVector<Scalar>&, int*);
   void subtract(GenVector<Scalar>&, int);
   void zeroAll();

   void putIn(Scalar *array, int position,  int num);
   void getFrom(Scalar *array, int position, int numdata);

   void getDataFrom(Scalar *array, int num);
   void addDataFrom(double *array, int num);

   void setn(int _n) {n = _n; }; 
   void setData(Scalar *v, int l) { len =l; d =v; }
   void setData(const GenVector<Scalar> &v1); 
   void insertData(Scalar *v);
   Scalar* getData() { return d; }
   void reset(int newlen, Scalar initialValue = 0.0);
   void resize(int newlen); // no-op if the sizes match, otherwise data is lost
   void conservativeResize(int newlen); // resizing with data preservation

   double squareNorm() const;
   double norm() const;
   double sqNorm() const { return ScalarTypes::norm((*this) * (*this)); }

   typedef int InfoType;
   typedef Scalar DataType;

//Extra operators used in non-linear elements----FL----2002----//

   void  vectorToTensor(Tensor_d1s0 &t) ;

//----------------------------------------------------------------------
  void setBlockValue(int k, int i, Scalar s) { d[k*n + i] = s; }
//  Scalar* getBlock(int k) { return d + n*k; }
  GenVector<Scalar>& getBlock(int k);
  void scaleBlock(int k, Scalar s) { for(int i=0; i<n; ++i) d[k*n+i] /= s; }
  void computeBlockNorms() {std::cerr << "GenVector::computeBlockNorms() called" << std::endl;}
  void printBlockNorms() {std::cerr << "GenVector::printBlockNorms() called" << std::endl;}
  double* getBlockNorms() { std::cerr << "GenVector::getBlockNorms() called" << std::endl; return 0;}
  void setNnzBlocks(int* bl) {};
  void printBlockDetails() {};
  int isnnz(int i) {return 1;}
};


template<class Scalar>
double norm(const GenVector<Scalar> &v) { return v.norm(); }

template<class Scalar>
inline Scalar &
GenVector<Scalar>::operator[] (int i) const { return d[i]; }

template<class Scalar>
inline
GenVector<Scalar>::GenVector(int l)           // Constructor
{ 
   myMemory = true;
   len = l;
   d   = new Scalar[len];
}

template<class Scalar>
inline
void
GenVector<Scalar>::initialize(int l)
{
   clean_up();
   myMemory = true;
   len = l;
   d = new Scalar[len];
}

template<class Scalar>
inline
GenVector<Scalar>::GenVector(const SingleInfo &inf)
{
  myMemory = true;
  len = inf.totLen();
  d   = new Scalar[len];
}

template<class Scalar>
inline
GenVector<Scalar>::GenVector(const BlockInfo &inf)
{
  myMemory = true;
  len = inf.numblocks * inf.blocklen;
  n = inf.blocklen;
  d   = new Scalar[len];
}

template<class Scalar>
inline
GenVector<Scalar>::GenVector(const GenVector<Scalar> &v, Scalar val)
{
 myMemory = true;
 len = v.len;
 d   = new Scalar[len];
 copy(val);
}

template<class Scalar>
inline
GenVector<Scalar>::~GenVector()
{
 if(d && myMemory) { delete [] d; d = 0; } //CBM
}

template<class Scalar>
inline
void
GenVector<Scalar>::clean_up()
{
 if(d && myMemory) { delete [] d; d=0; }
}

template<class Scalar>
GenVector<Scalar>
operator*(const Scalar c,const GenVector<Scalar> &v);

template<class Scalar>
class GenStackVector : public GenVector<Scalar> {
  public:
    GenStackVector() { this->len = 0; this->d = 0; }
    GenStackVector(Scalar *ptr, int length);
    GenStackVector(int length, Scalar *ptr);
    ~GenStackVector() { this->d=0; }
};


template<class Scalar>
inline
GenStackVector<Scalar>::GenStackVector(Scalar *ptr, int length)
{
  this->len = length;
  this->d   = ptr;
}

template<class Scalar>
inline
GenStackVector<Scalar>::GenStackVector(int length, Scalar *ptr)
{
  this->len = length;
  this->d   = ptr;
}

typedef GenVector<double> Vector;
typedef GenStackVector<double> StackVector;
typedef GenVector<DComplex> ComplexVector;
typedef GenStackVector<DComplex> ComplexStackVector;

#ifdef _TEMPLATE_FIX_
#include <Math.d/Vector.C>
#endif

#endif
