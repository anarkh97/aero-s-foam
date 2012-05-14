#ifndef _NONLINEARTRANSLATIONALSPRING_H_
#define _NONLINEARTRANSLATIONALSPRING_H_

#include <Element.d/Joint.d/DotConstraintType2.h>

// this element is a translational spring for large displacements and rotations, implemented using
// the penalized constraint method, in which the penalty parameter is the value of the spring stiffness.
// The constraint which is penalized in the scalar projection of vector d (from node A to node B)
// onto the unit vector c1 (either the x,y or z axis of the local frame attached to node A) minus 
// the constant value of same scalar projection in the undeformed configuration.
// i.e f(x) = d.dot(c1) - sp0

class NonlinearTranslationalSpring : public DotConstraintType2
{
    double sp0; // scalar projection of d onto c0 in the undeformed configuration

  public:
    NonlinearTranslationalSpring(int*, int);
    void setProp(StructProp *p, bool _myProp = false);
    void buildFrame(CoordSet&);
    void update(GeomState& gState, CoordSet& cs, double);
};

#endif
