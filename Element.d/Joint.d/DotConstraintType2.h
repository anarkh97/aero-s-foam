#ifndef _DOTCONSTRAINTTYPE2_H_
#define _DOTCONSTRAINTTYPE2_H_

#include <Element.d/MpcElement.d/ConstraintElement.h>

class DotConstraintType2 : public ConstraintElement
{
    EFrame *elemframe;
    int axis;
    double c0[3][3]; // initial frame (axes stored row-wise)
  public:
    DotConstraintType2(int*, int);
    void computeMPCs(CoordSet&);
    int getTopNumber();
    void updateLMPCs(GeomState& gState, CoordSet& cs);
    void getHessian(GeomState& gState, CoordSet&, FullSquareMatrix& H);
    void setFrame(EFrame *);
};

#endif
