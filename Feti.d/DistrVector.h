#ifndef _DISTRVECTOR_H_
#define _DISTRVECTOR_H_

#include <stdio.h>
#include <Driver.d/Communicator.h>
#include <Utils.d/MyComplex.h>

struct DistrInfo {

   int len;
   union {
     int numDom;
     int numLocSub;
   };
   union {
     int *domLen;
     int *subLen;
   };
   int *subOffset;
   // For parallel operations grouped by thread:
   int numLocThreads;
   int *threadOffset;
   int *threadLen;
   bool *masterFlag;
         
   FSCommunicator *com;
   
   DistrInfo(int i);
   DistrInfo() { initialize(); };
   DistrInfo(const DistrInfo &d);
   ~DistrInfo(); 
   void setMasterFlag();
   void setMasterFlag(bool *_masterFlag) { masterFlag = _masterFlag; }
   void computeOffsets();
   int totLen() { return len; }
   int *getMasterFlag(int i) { return 0; }
 private:
   void initialize();
};
/*
//------------------------------------------------------------------------------

template<class T, class Scalar, class IType = typename T::InfoType>
class Expr {

public:

  typedef IType InfoType;
  InfoType inf;
  T x;

  Expr(T v) : x(v), inf(x.info()) {}
  Expr(T v, IType i) : x(v), inf(i) { }

  Scalar operator[] (int i) const { return x[i]; }
  IType info() const { return inf; }

};
*/
#include <Math.d/Expr.h>
//------------------------------------------------------------------------------

template<class Scalar> class GenPartialDistrVector;

template<class Scalar>
class GenDistrVector {
  protected:
    bool myMemory;
    int len;		// entire length of the vector
    int numDom;		// number of domains
    Scalar *v;		// entire vector
    Scalar **subV;	// pointers to each domains sub-vector
    int *subVLen;	// length of each domains sub-vector
    int nT, *thLen;     // number of threads and lengths per thread
    Scalar **thV;       // each thread's subvector
    int *subVOffset, *thOffset;
    bool *masterFlag;
    bool infoFlag;
    const DistrInfo &inf;
    Scalar *partial;
  public:
    GenDistrVector() : inf(*(new DistrInfo)) { len = numDom = nT = 0; } 
    //GenDistrVector() { len = numDom = nT = 0; } // YYY was commented out and giving problem in initialization in Sfem
    GenDistrVector(const DistrInfo &dinfo);
    GenDistrVector(const GenDistrVector<Scalar> &v);
    GenDistrVector(const DistrInfo &dinfo, Scalar *, bool myMemory = true);
    virtual ~GenDistrVector();
    void initialize();
    void zero();
    void clean_up();
    int size() const { return len; }
    int num() const { return numDom; }
    Scalar &operator[](int i) { return v[i]; }
    Scalar operator[](int i) const { return v[i]; } 
    Scalar operator * (GenDistrVector &);
    Scalar operator ^ (GenDistrVector &);
    double norm();
    double infNorm();
    double sqNorm() { return ScalarTypes::norm((*this) * (*this)); }
    GenDistrVector &operator=(GenDistrVector<Scalar> &);
    GenDistrVector &operator=(Scalar c);
    GenDistrVector &operator*=(Scalar c);
    GenDistrVector &operator/=(Scalar c);
    GenDistrVector &operator+=(GenDistrVector<Scalar> &);
    GenDistrVector &operator-=(GenDistrVector<Scalar> &);
    GenDistrVector &operator/=(GenDistrVector<Scalar> &);
    GenDistrVector &linAdd(GenDistrVector<Scalar> &);
    GenDistrVector &linAdd(Scalar, GenDistrVector<Scalar> &);
    GenDistrVector &linAdd_inv(Scalar, GenDistrVector<Scalar> &);
    GenDistrVector &linAdd(Scalar, GenPartialDistrVector<Scalar> &);
    GenDistrVector &linAdd(Scalar, GenDistrVector<Scalar> &, Scalar, GenDistrVector<Scalar> &);
    GenDistrVector &linC(GenDistrVector<Scalar> &, Scalar);
    GenDistrVector &linC(GenDistrVector<Scalar> &, Scalar, GenDistrVector<Scalar> &);
    GenDistrVector &linC(Scalar, GenDistrVector<Scalar> &, Scalar, GenDistrVector<Scalar> &);
    GenDistrVector &swap(GenDistrVector<Scalar> &);
    template <class T>
      GenDistrVector &operator=(const Expr<T,Scalar> &);
    template <class T>
      GenDistrVector &operator+=(const Expr<T,Scalar> &);
    template <class T>
      GenDistrVector &operator-=(const Expr<T,Scalar> &);

  //  GenDistrVector operator+(GenDistrVector<Scalar> &);
  //  GenDistrVector operator-(GenDistrVector<Scalar> &);
    void updateBlock(int ii, Scalar c, GenDistrVector<Scalar> &) {cerr << "GenDistrVector::updateBlock not implemented" << endl;}
    void copyBlock(GenDistrVector<Scalar> &, int ii) {cerr << "GenDistrVector::copyBlock not implemented" << endl;}
    void addBlockSqr(int ii, Scalar c, GenDistrVector<Scalar> &);
    void computeSqrt();
    void computeRealz(int ii, Scalar c, GenDistrVector<Scalar> &) {cerr << "GenDistrVector::computeRealz not implemented" << endl;}
    void setn(int _n) {};      
    GenDistrVector<Scalar>&  getBlock(int iblock) { cerr << "GenDistrVector::getBlock not implemented" << endl; 
                                                    return *(new GenDistrVector<Scalar>()); }

    void negate();
    Scalar *data() const      { return v;          }
    Scalar *subData(int i)    { return subV[i];    }
    int subLen(int i)         { return subVLen[i]; }
    int numThreads()	      { return nT; }
    int threadLen(int i)      { return thLen[i];   }
    int threadOffset(int i)   { return thOffset[i]; }
    Scalar *threadData(int i) { return thV[i];     }
    bool *threadMasterFlag(int i) { return masterFlag +thOffset[i]; }
    bool *subMasterFlag(int i) { return masterFlag +subVOffset[i]; }
    Scalar ident();
    void print();
    void printNonZeroTerms();
    void printAll();
    void initRand();
    void doubleUp(int, Scalar, GenDistrVector<Scalar> * , GenDistrVector<Scalar> * , 
                  GenDistrVector<Scalar> *);
    void tripleUp(int, Scalar, GenDistrVector<Scalar> * , GenDistrVector<Scalar> * ,
                  GenDistrVector<Scalar> *, GenDistrVector<Scalar> * , GenDistrVector<Scalar> *);
    virtual Scalar getPartial(int iSub) { return 1.0; }

    Scalar sum() {
       Scalar x =0;
       for(int i = 0; i < len; ++i) x+= v[i];
       return x;
    }

   void scaleBlock(int k, Scalar s) { cerr << "Error : GenDistrVector::scaleBlock not implemented " << endl; } 

    typedef const DistrInfo &InfoType;
    const DistrInfo &info() const { return inf; }
};

template<class Scalar>
class GenStackDistVector : public GenDistrVector<Scalar> {
   public:
      GenStackDistVector(const DistrInfo &dinfo, Scalar *v)
               : GenDistrVector<Scalar>(dinfo, v, false) { }
      virtual ~GenStackDistVector() { }
};

template<class Scalar>
class GenPartialDistrVector : public GenDistrVector<Scalar> {
     // PJSA 1-22-07 distributed vector for which some subdomains have zero subvectors
     // can speed up some operations like * by skipping these subvectors
     //Scalar *partial;
   public:
     GenPartialDistrVector(const DistrInfo &dinfo)
        : GenDistrVector<Scalar>(dinfo) { this->partial = new Scalar[this->numDom]; for(int i=0; i<this->numDom; ++i) this->partial[i] = 1.0; }
     virtual ~GenPartialDistrVector() { delete [] this->partial; }

     Scalar operator * (GenDistrVector<Scalar> &);
     void computePartial();
     Scalar getPartial(int iSub) { return this->partial[iSub]; }
};

typedef GenDistrVector<double> DistrVector;
typedef GenDistrVector<DComplex> ComplexDistrVector;
typedef GenStackDistVector<double> StackDistVector;
typedef GenStackDistVector<DComplex> ComplexStackDistVector;

//------------------------------------------------------------------------------
/*
template<class A, class B>
class ProdRes {
  public:
    typedef B ResType;
};


template<>
class ProdRes<complex<double>, double> {
  public:
    typedef complex<double> ResType;
};
*/
//-----------------------------------------------------------------------------
template<class T1, class T2, class Scalar>
Scalar operator,(const Expr<T1,Scalar> &,const Expr<T2,Scalar> &);
//-----------------------------------------------------------------------------
template<class T2, class Scalar>
Scalar operator,(const GenDistrVector<Scalar> &,const Expr<T2,Scalar> &);
//-----------------------------------------------------------------------------
template<class T1,class Scalar>
Scalar operator,(const Expr<T1,Scalar> &,const GenDistrVector<Scalar> &);
//-----------------------------------------------------------------------------
template<class Scalar>
Scalar operator,(const GenDistrVector<Scalar> &,const GenDistrVector<Scalar> &);
//------------------------------------------------------------------------------
template<class Scalar>
double norm(const GenDistrVector<Scalar> &);
//------------------------------------------------------------------------------
template<class T1, class Scalar>
double norm(const Expr<T1,Scalar,const DistrInfo&> &);
//------------------------------------------------------------------------------
/*
template<class T1, class T2, class Scalar, class IType = typename T1::InfoType>
class Sum {
public:
  typedef IType InfoType;
private:
  T1 a;
  T2 b;
  InfoType len;

public:

  Sum(T1 aa, T2 bb, InfoType l) : a(aa), b(bb), len(l) { }

  Scalar operator[](int i) const { return a[i]+b[i]; }
  InfoType info() const { return len; }

};
//------------------------------------------------------------------------------

template<class T1, class T2, class Scalar>
inline
Expr<Sum<T1, T2, Scalar>, Scalar>
operator+(const Expr<T1, Scalar> &x1, const Expr<T2, Scalar> &x2)
{

  return Expr<Sum<T1, T2, Scalar>, Scalar>
    ( Sum<T1, T2, Scalar>(x1.x, x2.x, x1.info()) );

}
*/
//------------------------------------------------------------------------------

template<class Scalar>
inline
Expr<
 Sum<const GenDistrVector<Scalar> &, const GenDistrVector<Scalar> &, Scalar,
    typename GenDistrVector<Scalar>::InfoType>
    , Scalar>
operator+(const GenDistrVector<Scalar> &v1, const GenDistrVector<Scalar> &v2)
{

  return Expr<
    Sum<const GenDistrVector<Scalar> &, const GenDistrVector<Scalar> &, Scalar,
     typename GenDistrVector<Scalar>::InfoType> , Scalar >
    ( Sum<const GenDistrVector<Scalar> &, const GenDistrVector<Scalar> &,
      Scalar, typename GenDistrVector<Scalar>::InfoType>(v1, v2, v1.info()) );

}

//------------------------------------------------------------------------------

template<class T, class Scalar>
inline
Expr<Sum<T, Scalar *, Scalar>, Scalar>
operator+(const Expr<T, Scalar> &x, const GenDistrVector<Scalar> &v)
{

  return Expr<Sum<T, Scalar *, Scalar>, Scalar>
    ( Sum<T, Scalar *, Scalar>(x.x, v.data(), v.info()) );

}

//------------------------------------------------------------------------------

template<class T, class Scalar>
inline
Expr<Sum<const GenDistrVector<Scalar> &, T, Scalar, 
   typename GenDistrVector<Scalar>::InfoType>, Scalar>
operator+(const GenDistrVector<Scalar> &v, const Expr<T, Scalar> &x)
{

  return Expr<Sum<const GenDistrVector<Scalar> &, T, Scalar, 
   typename GenDistrVector<Scalar>::InfoType>
   , Scalar>
    ( Sum<const GenDistrVector<Scalar> &, T, Scalar, 
   typename GenDistrVector<Scalar>::InfoType>(v.data(), x.x, v.info()) );

}

//------------------------------------------------------------------------------
/*
template<class T1, class T2, class Scalar, class IType = typename T1::InfoType>
class Diff {
public:
  typedef IType InfoType;
private:
  T1 a;
  T2 b;
  InfoType len;

public:

  Diff(T1 aa, T2 bb, InfoType l) : a(aa), b(bb), len(l) { }

  Scalar operator[](int i) const { return a[i]+b[i]; }
  InfoType info() const { return len; }

};
//-----------------------------------------------------------------------------

//------------------------------------------------------------------------------

template<class T1, class T2, class Scalar>
inline
Expr<Diff<T1, T2, Scalar>, Scalar>
operator-(const Expr<T1, Scalar> &x1, const Expr<T2, Scalar> &x2)
{

  return Expr<Diff<T1, T2, Scalar>, Scalar>
    ( Diff<T1, T2, Scalar>(x1.x, x2.x, x1.info()) );

}
*/
//------------------------------------------------------------------------------

template<class Scalar>
inline
Expr<
 Diff<const GenDistrVector<Scalar> &, const GenDistrVector<Scalar> &, Scalar,
    typename GenDistrVector<Scalar>::InfoType>
    , Scalar>
operator-(const GenDistrVector<Scalar> &v1, const GenDistrVector<Scalar> &v2)
{

  return Expr<
    Diff<const GenDistrVector<Scalar> &, const GenDistrVector<Scalar> &, Scalar,
     typename GenDistrVector<Scalar>::InfoType> , Scalar >
    ( Diff<const GenDistrVector<Scalar> &, const GenDistrVector<Scalar> &,
      Scalar, typename GenDistrVector<Scalar>::InfoType>(v1, v2, v1.info()) );

}

//------------------------------------------------------------------------------

template<class T, class Scalar>
inline
Expr<Diff<T, const GenDistrVector<Scalar> &, Scalar>, Scalar>
operator-(const Expr<T, Scalar> &x, const GenDistrVector<Scalar> &v)
{

  return Expr<Diff<T, const GenDistrVector<Scalar> &, Scalar>, Scalar>
    ( Diff<T, const GenDistrVector<Scalar> &, Scalar>(x.x, v, v.info()) );

}

//------------------------------------------------------------------------------

template<class T, class Scalar>
inline
Expr<Diff<const GenDistrVector<Scalar> &, T, Scalar, 
   typename GenDistrVector<Scalar>::InfoType>, Scalar>
operator-(const GenDistrVector<Scalar> &v, const Expr<T, Scalar> &x)
{

  return Expr<Diff<const GenDistrVector<Scalar> &, T, Scalar, 
   typename GenDistrVector<Scalar>::InfoType>
   , Scalar>
    ( Diff<const GenDistrVector<Scalar> &, T, Scalar, 
   typename GenDistrVector<Scalar>::InfoType>(v.data(), x.x, v.info()) );

}

//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
/*
template<class T, class Scalar, class Res, class IType = typename T::InfoType>
class OuterProd {

public:
  typedef IType InfoType;
private:
  Scalar y;
  T a;
  InfoType len;

public:

  OuterProd(Scalar yy, T aa, InfoType l) : y(yy), a(aa), len(l) { }

  Res operator[](int i) const { return y*a[i]; }
  InfoType info() const { return len; }

};

//------------------------------------------------------------------------------

template<class T, class S2>
inline
Expr<OuterProd<T, double, typename ProdRes<double, S2>::ResType>,
    typename ProdRes<double,S2>::ResType > operator*(double y, const Expr<T, S2> &x)
{

  return Expr<OuterProd<T, double, typename ProdRes<double,S2>::ResType>,
         typename ProdRes<double,S2>::ResType>
    ( OuterProd<T, double, typename ProdRes<double,S2>::ResType>(y, x.x, x.info()) );
}

//------------------------------------------------------------------------------

template<class T, class S2>
inline
Expr<OuterProd<T, complex<double>, typename ProdRes<complex<double>, S2>::ResType>,
    typename ProdRes<complex<double>,S2>::ResType > operator*(complex<double> y, const Expr<T, S2> &x)
{

  return Expr<OuterProd<T, complex<double>, typename ProdRes<complex<double>,S2>::ResType>,
         typename ProdRes<complex<double>,S2>::ResType>
    ( OuterProd<T, complex<double>, typename ProdRes<complex<double>,S2>::ResType>(y, x.x, x.info()) );
}
*/
//------------------------------------------------------------------------------
template<class Scalar, class Res>
inline
Expr<OuterProd<const GenDistrVector<Res> &, Scalar, 
       typename ProdRes<Scalar,Res>::ResType, 
       typename GenDistrVector<Res>::InfoType>,
     typename ProdRes<Scalar,Res>::ResType> operator*(Scalar y, const
     GenDistrVector<Res> &v)
{

  return Expr<
     OuterProd<const GenDistrVector<Res> &, Scalar, 
       typename ProdRes<Scalar,Res>::ResType, 
       typename GenDistrVector<Res>::InfoType>,
         typename ProdRes<Scalar,Res>::ResType>
    ( OuterProd<const GenDistrVector<Res> &, Scalar, 
       typename ProdRes<Scalar,Res>::ResType,
       typename GenDistrVector<Res>::InfoType>(y,
    v, v.info()) );

}

class CompIsL {
  public:
    template <class A, class B>
    static bool apply(const A &a, const B&b) { return a < b; }
};

class CompIsLE {
  public:
    template <class A, class B>
    static bool apply(const A &a, const B&b) { return a <= b; }
};

class CompIsG {
  public:
    template <class A, class B>
    static bool apply(const A &a, const B&b) { return a > b; }
};

class CompIsGE {
  public:
    template <class A, class B>
    static bool apply(const A &a, const B&b) { return a >= b; }
};

class CompIsE {
  public:
    template <class A, class B>
    static bool apply(const A &a, const B&b) { return a == b; }
};

class BoolAnd {
  public:
    template <class A, class B>
    static bool apply(const A &a, const B&b) { return a && b; }
};

class BoolOr {
  public:
    template <class A, class B>
    static bool apply(const A &a, const B&b) { return a || b; }
};

/*
template<class T1, class T2, class CompOp, class IType = typename T1::InfoType>
class CmpElem {
public:
  typedef IType InfoType;
private:
  T1 a;
  T2 b;
  InfoType len;

public:

  CmpElem(T1 aa, T2 bb, InfoType l) : a(aa), b(bb), len(l) { }

  bool operator[](int i) const { return CompOp::apply(a[i],b[i]); }
  InfoType info() const { return len; }

};
*/

// macro calls to create the operators
binaryOpDec( == , CompIsE , , GenDistrVector<Scalar> )
/*
template<class Scalar>
bool operator==(const  GenDistrVector<Scalar> &, const  GenDistrVector<Scalar>
&);*//*
template<class Scalar >
 inline Expr< 
   BinOp<const GenDistrVector<Scalar> &, const GenDistrVector<Scalar> &,
      Scalar, CompIsE,  typename GenDistrVector<Scalar>::InfoType>,
       Scalar > 
operator==(const GenDistrVector<Scalar> &v1, const GenDistrVector<Scalar> &v2)
 { return Expr<
 BinOp<const GenDistrVector<Scalar> &, const GenDistrVector<Scalar> &, 
 Scalar, CompIsE, typename GenDistrVector<Scalar>::InfoType>, 
 Scalar > ( BinOp<const GenDistrVector<Scalar> &,
  const GenDistrVector<Scalar> &, Scalar, CompIsE, 
  typename GenDistrVector<Scalar>::InfoType>
  (v1,v2)); }*/

//-----------------------------------------------------------------------------

template<class T, class IType = typename T::InfoType>
class BoolNot {
public:
  typedef IType InfoType;
private:
  T a;
  InfoType len;
public:
  BoolNot(T t, IType it) : a(t), len(it) {}
  bool operator[](int i) const { return !a[i]; }
  InfoType info() const { return len; }

};
#ifdef _TEMPLATE_FIX_
#include <Feti.d/DistrVector.C>
#endif

#endif
