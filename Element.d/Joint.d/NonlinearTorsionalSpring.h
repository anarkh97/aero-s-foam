#ifndef _NONLINEARTORSIONALSPRING_H_
#define _NONLINEARTORSIONALSPRING_H_

#include <Element.d/MpcElement.d/MpcElement.h>

// this element is a torsional spring for large displacements and rotations, implemented using
// the penalized constraint method, in which the penalty parameter is the value of the spring stiffness.
// The constraint which is penalized in the relative rotation of node A w.r.t node B about a specified axis.
// i.e f(x) = arccos(c1.dot(c2))
// To handle arbitriarily large rotations and also avoid singularities dynamic reparameterization and 
// history variables are employed

class NonlinearTorsionalSpring : public MpcElement
{
    double (*c0)[3]; // initial frame (axes stored row-wise)
    int m_axis1, m_axis2;
    bool covariant_derivatives;

    double offset, offset2;
    int quadrant;

  public:
    NonlinearTorsionalSpring(int*, int, int);
    ~NonlinearTorsionalSpring();
    void setProp(StructProp *p, bool _myProp = false);
    void setFrame(EFrame *);
    void buildFrame(CoordSet&);
    void update(GeomState& gState, CoordSet& cs, double t);
    void getHessian(GeomState& gState, CoordSet&, FullSquareMatrix& H, double);
    int numStates();
    void initStates(double *);
    void updateStates(GeomState *refState, GeomState &curState, CoordSet &C0);
};

#endif
