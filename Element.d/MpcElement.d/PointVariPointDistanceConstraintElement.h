#ifndef _POINTVARIPOINTDISTANCECONSTRAINTELEMENT_H_
#define _POINTVARIPOINTDISTANCECONSTRAINTELEMENT_H_

#include <Element.d/Function.d/Constraint.d/PointVariPointDistanceConstraintFunction.h>
#include <Element.d/MpcElement.d/ConstraintFunctionElement.h>

class PointVariPointDistanceConstraintElement : public ConstraintFunctionElement<Simo::PointVariPointDistanceConstraintFunction>
{
  public:
    PointVariPointDistanceConstraintElement(int* _nn); 

  protected:
    void getConstants(CoordSet& cs, Eigen::Array<double,11,1>& sconst, Eigen::Array<int,1,1>& iconst, GeomState* = NULL);
};

#endif
