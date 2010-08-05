#ifndef _DOTCONSTRAINTTYPE1A_H_
#define _DOTCONSTRAINTTYPE1A_H_

#include <Element.d/MpcElement.d/MpcElement.h>

class DotConstraintType1a : public MpcElement
{
    EFrame *elemframe;
    int axis1, axis2;
    double c0[3][3]; // initial frame (axes stored row-wise)
    double t; // time
  public:
    DotConstraintType1a(int*, int, int);
    void buildFrame(CoordSet&);
    int getTopNumber();
    void update(GeomState& gState, CoordSet& cs, double t);
    void getHessian(GeomState& gState, CoordSet&, FullSquareMatrix& H);
    void setFrame(EFrame *);
};

#endif
