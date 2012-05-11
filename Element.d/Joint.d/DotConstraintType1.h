#ifndef _DOTCONSTRAINTTYPE1_H_
#define _DOTCONSTRAINTTYPE1_H_

#include <Element.d/MpcElement.d/MpcElement.h>

// Ref: Rigid Body Dynamics of Mechanisms Vol 1, Hubert Hahn, section 5.2.1.4
// Constraint Building Block type BB4, also known as a rotation blocker constraint
// one constrained rotational DOF

class DotConstraintType1 : public MpcElement
{
    double (*c0)[3]; // initial frame (axes stored row-wise)
    int axis1, axis2;
    bool covariant_derivatives;

  public:
    DotConstraintType1(int*, int, int);
    ~DotConstraintType1();
    void setFrame(EFrame *);
    void buildFrame(CoordSet&);
    void update(GeomState& gState, CoordSet& cs, double);
    void getHessian(GeomState& gState, CoordSet&, FullSquareMatrix& H);
};

#endif
