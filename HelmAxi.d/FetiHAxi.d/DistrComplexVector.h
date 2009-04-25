#ifndef _DISTRCOMPLEXVECTOR_H_
#define _DISTRCOMPLEXVECTOR_H_


#include <Utils.d/MyComplex.h>


class DistrInfo;


class DistrComplexVector {

    int len;           // entire length of the vector
    int numDom;        // number of domains
    DComplex *v;       // entire vector
    DComplex **subV;   // pointers to each subdomains sub-vector
    int *subVLen;      // length of each domains subvector
    int nT, *thLen;    // number of threads and lengths per thread
    DComplex **thV;    // each thread's subvector

  public:

    DistrComplexVector();
    DistrComplexVector(const DistrInfo &);

    void zero();

    DComplex operator * (DistrComplexVector&);
    DComplex operator ^ (DistrComplexVector&);
    DComplex norm();

    DistrComplexVector &operator=(DistrComplexVector &);
    DistrComplexVector &operator*=(DComplex c);
    DistrComplexVector &operator^=(DComplex c);
    DistrComplexVector &operator+=(DistrComplexVector &);
    DistrComplexVector &operator-=(DistrComplexVector &);
    DistrComplexVector &linAdd(DComplex, DistrComplexVector&);
    DistrComplexVector &linAdd(DComplex, DistrComplexVector&, 
                               DComplex, DistrComplexVector&);
    DistrComplexVector &linC(DistrComplexVector & , DComplex , 
                             DistrComplexVector & );
    DistrComplexVector &linC(DComplex, DistrComplexVector & , 
                             DComplex , DistrComplexVector & );

    DistrComplexVector &copy(DistrComplexVector &);

    void negate();
    DComplex *data()            { return v; }
    DComplex *subData(int i)    { return subV[i]; }
    int subLen(int i)           { return subVLen[i]; }
    int threadLen(int i)        { return thLen[i]; }
    DComplex *threadData(int i) { return thV[i]; }
    DComplex ident();
    void print();
    void printAll();

    friend class VecOpComplex;
};

#endif
