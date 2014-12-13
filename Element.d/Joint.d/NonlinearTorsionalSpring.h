#ifndef _NONLINEARTORSIONALSPRING_H_
#define _NONLINEARTORSIONALSPRING_H_

#include <Element.d/MpcElement.d/AngleType1ConstraintElement.h>

// this element is a torsional spring for large displacements and rotations, implemented using
// the penalized constraint method, in which the penalty parameter is the value of the spring stiffness.
// The constraint which is penalized in the relative rotation of node A w.r.t node B about a specified axis.
// i.e f(x) = arccos(c1.dot(c2))
// To handle arbitriarily large rotations and also avoid singularities dynamic reparameterization and 
// history variables are employed

class NonlinearTorsionalSpring : public AngleType1ConstraintElement
{
    int m_axis1, m_axis2;
    double offset2;
    int quadrant;

  public:
    NonlinearTorsionalSpring(int*, int, int);
    void setProp(StructProp *p, bool _myProp = false);
    void update(GeomState *refState, GeomState& gState, CoordSet& cs, double t);

    int numStates();
    void initStates(double *);
    void updateStates(GeomState *refState, GeomState &curState, CoordSet &C0, double dt = 0);

    bool isSpring() { return true; }
    bool hasRot() { return true; }
};

#endif
