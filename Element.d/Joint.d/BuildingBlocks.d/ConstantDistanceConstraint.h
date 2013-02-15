#ifndef _CONSTANTDISTANCECONSTRAINT_H_
#define _CONSTANTDISTANCECONSTRAINT_H_

#include <Element.d/MpcElement.d/DistanceConstraintElement.h>

// Ref: Rigid Body Dynamics of Mechanisms Vol 1, Hubert Hahn, sections 5.2.1.5 and 5.2.2.2
// Constraint Building Block type BB5, also known as massless spherical-spherical link
// one constrained translational DOF

class ConstantDistanceConstraint : public DistanceConstraintElement
{
  public:
    ConstantDistanceConstraint(int*);
    void buildFrame(CoordSet&);
    int getTopNumber();
};

#endif