#ifndef _FLEXIBLETWONODETRUSS_H_
#define _FLEXIBLETWONODETRUSS_H_

#include <Element.d/SuperElement.h>

class FlexibleTwoNodeTruss : public SuperElement
{
    double l0;

  public:
    FlexibleTwoNodeTruss(int*);
    int getTopNumber() { return 101; }
    bool isSafe() { return false; }
    PrioInfo examine(int sub, MultiFront*);

    void buildFrame(CoordSet&);
    void setProp(StructProp *p, bool myProp);
    int getMassType() { return 2; }
    FullSquareMatrix massMatrix(CoordSet &cs, double *mel, int cmflg = 1);
    double getMass(CoordSet& cs);
};

#endif
