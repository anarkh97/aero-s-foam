#ifndef _DOTCONSTRAINTTYPE2A_H_
#define _DOTCONSTRAINTTYPE2A_H_

#include <Element.d/MpcElement.d/DotType2ConstraintElement.h>

class DotConstraintType2a : public DotType2ConstraintElement
{
  public:
    DotConstraintType2a(int*, int);
    void buildFrame(CoordSet& cs);
    void update(GeomState& gState, CoordSet& cs, double t);
    double getVelocityConstraintRhs(GeomState *gState, CoordSet& cs, double t);
    double getAccelerationConstraintRhs(GeomState *gState, CoordSet& cs, double t);
    void computePressureForce(CoordSet&, Vector& elPressureForce,
                              GeomState *gs = 0, int cflg = 0, double t = 0);
};

#endif
