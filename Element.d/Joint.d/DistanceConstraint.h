#ifndef _DISTANCECONSTRAINT_H_
#define _DISTANCECONSTRAINT_H_

#include <Element.d/MpcElement.d/ConstraintElement.h>

class DistanceConstraint : public ConstraintElement
{
    EFrame *elemframe;
    double l0;       // initial length
  public:
    DistanceConstraint(int*);
    void computeMPCs(CoordSet&);
    int getTopNumber();
    void updateLMPCs(GeomState& gState, CoordSet& cs);
    void getHessian(GeomState& gState, CoordSet&, FullSquareMatrix& H);
};

#endif
