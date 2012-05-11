#ifndef _TORSIONALSPRINGTYPE1_H_
#define _TORSIONALSPRINGTYPE1_H_

#include <Element.d/MpcElement.d/MpcElement.h>

class TorsionalSpringType1 : public MpcElement
{
    double (*c0)[3]; // initial frame (axes stored row-wise)
    int m_axis1, m_axis2;
    bool covariant_derivatives;

    double offset, offset2;
    int quadrant;

  public:
    TorsionalSpringType1(int*, int, int);
    ~TorsionalSpringType1();
    void setProp(StructProp *p, bool _myProp = false);
    void setFrame(EFrame *);
    void buildFrame(CoordSet&);
    void update(GeomState& gState, CoordSet& cs, double t);
    void getHessian(GeomState& gState, CoordSet&, FullSquareMatrix& H);
    int numStates();
    void updateStates(GeomState *refState, GeomState &curState, CoordSet &C0);
};

#endif
