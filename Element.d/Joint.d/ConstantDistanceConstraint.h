#ifndef _CONSTANTDISTANCECONSTRAINT_H_
#define _CONSTANTDISTANCECONSTRAINT_H_

#include <Element.d/MpcElement.d/MpcElement.h>

class ConstantDistanceConstraint : public MpcElement
{
    double l0;       // initial length
  public:
    ConstantDistanceConstraint(int*);
    void buildFrame(CoordSet&);
    int getTopNumber();
    void update(GeomState& gState, CoordSet& cs, double);
    void getHessian(GeomState& gState, CoordSet&, FullSquareMatrix& H);
};

#endif
