#ifndef _SPHERICALJOINT_H_
#define _SPHERICALJOINT_H_

#include <Element.d/SuperElement.h>

// Ref: Rigid Body Dynamics of Mechanisms Vol 1, Hubert Hahn, sections 5.2.1.1 and 5.2.2.1
// Constraint Building Block type BB1, also known as a common point constraint
// constrains three translational dofs

class SphericalJoint : public SuperElement
{
  public:
    SphericalJoint(int*);
    int getTopNumber();
    PrioInfo examine(int sub, MultiFront*);
};

#endif
