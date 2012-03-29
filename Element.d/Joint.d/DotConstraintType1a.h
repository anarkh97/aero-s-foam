#ifndef _DOTCONSTRAINTTYPE1A_H_
#define _DOTCONSTRAINTTYPE1A_H_

#include <Element.d/MpcElement.d/MpcElement.h>

class DotConstraintType1a : public MpcElement
{
    double (*c0)[3]; // initial frame (axes stored row-wise)
    int axis1, axis2;
    int axis1_copy;
    double t_reparam, offset;
    bool instantaneous_partials;

  public:
    DotConstraintType1a(int*, int, int);
    ~DotConstraintType1a();
    void setFrame(EFrame *);
    void buildFrame(CoordSet&);
    void update(GeomState& gState, CoordSet& cs, double t);
    void getHessian(GeomState& gState, CoordSet&, FullSquareMatrix& H);
    double getVelocityConstraintRhs(GeomState *gState, CoordSet& cs, double t);
    double getAccelerationConstraintRhs(GeomState *gState, CoordSet& cs, double t);
    void computePressureForce(CoordSet&, Vector& elPressureForce,
                              GeomState *gs = 0, int cflg = 0, double t = 0);
};

#endif
