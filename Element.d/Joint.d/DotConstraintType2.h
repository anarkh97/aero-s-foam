#ifndef _DOTCONSTRAINTTYPE2_H_
#define _DOTCONSTRAINTTYPE2_H_

#include <Element.d/MpcElement.d/MpcElement.h>

class DotConstraintType2 : public MpcElement
{
    EFrame *elemframe;
    int axis;
    double c0[3][3]; // initial frame (axes stored row-wise)
  public:
    DotConstraintType2(int*, int);
    void buildFrame(CoordSet&);
    int getTopNumber();
    void update(GeomState& gState, CoordSet& cs);
    void getHessian(GeomState& gState, CoordSet&, FullSquareMatrix& H);
    void setFrame(EFrame *);
};

#endif
