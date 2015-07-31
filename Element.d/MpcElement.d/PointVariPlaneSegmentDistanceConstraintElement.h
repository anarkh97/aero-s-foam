#ifndef _POINTVARIPLANESEGMENTDISTANCECONSTRAINTELEMENT_H_
#define _POINTVARIPLANESEGMENTDISTANCECONSTRAINTELEMENT_H_

#include <Element.d/Function.d/Constraint.d/PointVariPlaneSegmentDistanceConstraintFunction.h>
#include <Element.d/MpcElement.d/ConstraintFunctionElement.h>

class PointVariPlaneSegmentDistanceConstraintElement : public ConstraintFunctionElement<Simo::PointVariPlaneSegmentDistanceConstraintFunction>
{
  public:
    PointVariPlaneSegmentDistanceConstraintElement(int* _nn); 
  protected:
    void getConstants(CoordSet& cs, Eigen::Array<double,17,1>& sconst, Eigen::Array<int,1,1>& iconst, GeomState* = NULL);
};

#endif
