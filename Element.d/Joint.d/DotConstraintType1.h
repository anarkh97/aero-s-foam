#ifndef _DOTCONSTRAINTTYPE1_H_
#define _DOTCONSTRAINTTYPE1_H_

#include <Element.d/MpcElement.d/ConstraintElement.h>

class DotConstraintType1 : public ConstraintElement
{
    EFrame *elemframe;
    int axis1, axis2;
    double c0[3][3]; // initial frame (axes stored row-wise)
  public:
    DotConstraintType1(int*, int, int);
    void computeMPCs(CoordSet&);
    int getTopNumber();
    void updateLMPCs(GeomState& gState, CoordSet& cs);
    void getHessian(GeomState& gState, CoordSet&, FullSquareMatrix& H);
    void setFrame(EFrame *);
};

#endif
