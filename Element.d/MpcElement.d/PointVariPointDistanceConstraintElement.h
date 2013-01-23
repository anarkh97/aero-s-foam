#ifndef _POINTVARIPOINTDISTANCECFE_H_
#define _POINTVARIPOINTDISTANCECFE_H_

#include <Element.d/MpcElement.d/ConstraintFunction.d/PointVariPointDistanceConstraintFunction.h>
#include <Element.d/MpcElement.d/ConstraintFunctionElement.h>

class PointVariPointDistanceConstraintElement : public ConstraintFunctionElement<PointVariPointDistanceConstraintFunction>
{
  public:
    PointVariPointDistanceConstraintElement(int* _nn); 

  protected:
    void getConstants(CoordSet& cs, Eigen::Array<double,11,1>& sconst, Eigen::Array<int,1,1>& iconst, GeomState* = NULL);
};

#endif
