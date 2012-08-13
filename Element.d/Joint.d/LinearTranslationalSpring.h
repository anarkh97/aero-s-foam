#ifndef _LINEARTRANSLATIONALSPRING_H_
#define _LINEARTRANSLATIONALSPRING_H_

#include <Element.d/Joint.d/ConstantDistanceConstraint.h>

// this element is a translational spring for small displacements and rotations, implemented using
// the penalized constraint method, in which the penalty parameter is the value of the spring stiffness.
// The constraint which is penalized in the distance l from node A to node B minus
// the constant value of same distance in the undeformed configuration l0.
// i.e f(x) = l - l0

class LinearTranslationalSpring : public ConstantDistanceConstraint
{
  public:
    LinearTranslationalSpring(int*);
    void setProp(StructProp *p, bool _myProp = false);

    bool isSpring() { return true; }
};

#endif