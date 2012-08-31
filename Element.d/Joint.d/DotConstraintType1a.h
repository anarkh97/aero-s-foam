#ifndef _DOTCONSTRAINTTYPE1A_H_
#define _DOTCONSTRAINTTYPE1A_H_

#include <Element.d/MpcElement.d/DotType1ConstraintElement.h>

class DotConstraintType1a : public DotType1ConstraintElement
{
    int axis1_copy;
    double t_reparam, offset;

  public:
    DotConstraintType1a(int*, int, int, int=2);
    void buildFrame(CoordSet& cs);
    void update(GeomState& gState, CoordSet& cs, double t);
    double getVelocityConstraintRhs(GeomState *gState, CoordSet& cs, double t);
    double getAccelerationConstraintRhs(GeomState *gState, CoordSet& cs, double t);
    void computePressureForce(CoordSet&, Vector& elPressureForce,
                              GeomState *gs = 0, int cflg = 0, double t = 0);
};

#endif
