#ifndef _DOTCONSTRAINTTYPE1_H_
#define _DOTCONSTRAINTTYPE1_H_

#include <Element.d/MpcElement.d/MpcElement.h>

class DotConstraintType1 : public MpcElement
{
    EFrame *elemframe;
    int axis1, axis2;
    double c0[3][3]; // initial frame (axes stored row-wise)
  public:
    DotConstraintType1(int*, int, int);
    void buildFrame(CoordSet&);
    int getTopNumber();
    void update(GeomState& gState, CoordSet& cs, double);
    void getHessian(GeomState& gState, CoordSet&, FullSquareMatrix& H);
    void setFrame(EFrame *);
};

#endif
