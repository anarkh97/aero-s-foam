#ifndef _DOTCONSTRAINTTYPE2_H_
#define _DOTCONSTRAINTTYPE2_H_

#include <Element.d/MpcElement.d/MpcElement.h>

class DotConstraintType2 : public MpcElement
{
    double (*c0)[3]; // initial frame (axes stored row-wise)
    int axis;

  public:
    DotConstraintType2(int*, int);
    ~DotConstraintType2();
    void setFrame(EFrame *);
    void buildFrame(CoordSet&);
    void update(GeomState& gState, CoordSet& cs, double);
    void getHessian(GeomState& gState, CoordSet&, FullSquareMatrix& H);
};

#endif
