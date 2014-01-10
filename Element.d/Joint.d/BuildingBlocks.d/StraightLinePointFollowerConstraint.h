#ifndef _STRAIGHTLINEPOINTFOLLOWERCONSTRAINT_H_
#define _STRAIGHTLINEPOINTFOLLOWERCONSTRAINT_H_

#include <Element.d/SuperElement.h>

// Ref: Rigid Body Dynamics of Mechanisms Vol 1, Hubert Hahn, section 5.2.1.3
// Constraint Building Block type BB3, also known as a straight line point follower constraint
// two constrained translational dofs

class StraightLinePointFollowerConstraint : public SuperElement
{
  public:
    StraightLinePointFollowerConstraint(int*, int=0);
    int getTopNumber();
    bool hasRot() { return true; }
};

#endif
