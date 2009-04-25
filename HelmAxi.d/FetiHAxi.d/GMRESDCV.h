#ifndef _GMRES_DCV_H_
#define _GMRES_DCV_H_

#include <Utils.d/MyComplex.h>


class DistrInfo;
class DistrComplexVector;


class GMRESDCV {

    DistrInfo distrInfo;
    int numV;
    int maxV;
    DistrComplexVector **allV;
    DistrComplexVector **allW;
    DComplex *givensC;
    DComplex *givensS;
    DComplex *g;
    DComplex *y;
    DComplex *matrixH;
    double beta;

    void generateRotation(DComplex a, DComplex b, DComplex &cs, DComplex &ss);
    void applyRotation(DComplex &a, DComplex &b, DComplex cs, DComplex ss);

 public:

    GMRESDCV(int maxSize, DistrInfo & dI);
    ~GMRESDCV();
    void reInit();
    void init(DistrComplexVector & r0);
    void orthoAdd(DistrComplexVector & Fv, DistrComplexVector & v);
    double update();
    void solution(DistrComplexVector &u);
    void dumpH();
};

#endif
