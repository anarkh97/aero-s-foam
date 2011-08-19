#ifndef _POINTLINEDISTANCECONSTRAINT_H_
#define _POINTLINEDISTANCECONSTRAINT_H_

#include <Element.d/MpcElement.d/MpcElement.h>

class PointLineDistanceConstraint : public MpcElement
{
    double x1[3], x2[3]; // coordinates of the 2 points defining the line
    double denominator;
    double d0;           // initial distance from point to line
  public:
    PointLineDistanceConstraint(int*);
    void setFrame(EFrame *);
    void buildFrame(CoordSet&);
    int getTopNumber();
    void update(GeomState& gState, CoordSet& cs, double);
    void getHessian(GeomState& gState, CoordSet&, FullSquareMatrix& H);
};

#endif
