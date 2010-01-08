#ifndef _CONSTANTDISTANCECONSTRAINT_H_
#define _CONSTANTDISTANCECONSTRAINT_H_

#include <Element.d/MpcElement.d/ConstraintElement.h>

class ConstantDistanceConstraint : public ConstraintElement
{
    EFrame *elemframe;
    double l0;       // initial length
  public:
    ConstantDistanceConstraint(int*);
    void computeMPCs(CoordSet&);
    int getTopNumber();
    void update(GeomState& gState, CoordSet& cs);
    void getHessian(GeomState& gState, CoordSet&, FullSquareMatrix& H);
};

#endif
