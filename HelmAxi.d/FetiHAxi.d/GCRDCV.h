#ifndef _GCR_DCV_H_
#define _GCR_DCV_H_


#include <Utils.d/MyComplex.h>


class DistrInfo;
class DistrComplexVector;


class GCRDCV {

    DistrInfo distrInfo;
    int numP;
    int maxP;
    DComplex *allFPFP;
    DistrComplexVector **allP; 
    DistrComplexVector **allFP; 

 public:

    GCRDCV(int maxSize, DistrInfo & dI);
    ~GCRDCV();
    void orthoAdd(DistrComplexVector &, DistrComplexVector &, DComplex);
    void orthogonalize(DistrComplexVector &, DistrComplexVector &, 
                       DistrComplexVector &, DistrComplexVector &);
    int predict(DistrComplexVector &, DistrComplexVector &);
    int update(DistrComplexVector &, DistrComplexVector &, 
               DistrComplexVector &);
    int numDir() { return numP; }

};

#endif
