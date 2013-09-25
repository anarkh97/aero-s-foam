#ifndef _POINTPOINTDISTANCECONSTRAINTELEMENT_H_
#define _POINTPOINTDISTANCECONSTRAINTELEMENT_H_

#include <Element.d/Function.d/Constraint.d/PointPointDistanceConstraintFunction.h>
#include <Element.d/MpcElement.d/ConstraintFunctionElement.h>

class PointPointDistanceConstraintElement : public ConstraintFunctionElement<Simo::PointPointDistanceConstraintFunction>
{
    double x1[3]; // coordinates of the fixed point

  public:
    PointPointDistanceConstraintElement(int* _nn); 
    void setFrame(EFrame *);

  protected:
    void getConstants(CoordSet& cs, Eigen::Array<double,11,1>& sconst, Eigen::Array<int,1,1>& iconst, GeomState* = NULL);
};

#endif
