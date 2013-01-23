#ifndef _POINTVARILINEDISTANCECFE_H_
#define _POINTVARILINEDISTANCECFE_H_

#include <Element.d/MpcElement.d/ConstraintFunction.d/PointVariLineDistanceConstraintFunction.h>
#include <Element.d/MpcElement.d/ConstraintFunctionElement.h>

class PointVariLineDistanceConstraintElement : public ConstraintFunctionElement<PointVariLineDistanceConstraintFunction>
{
  public:
    PointVariLineDistanceConstraintElement(int* _nn); 

  protected:
    void getConstants(CoordSet& cs, Eigen::Array<double,14,1>& sconst, Eigen::Array<int,1,1>& iconst, GeomState* = NULL);
};

#endif
