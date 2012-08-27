#ifndef _POINTVARIPLANEDISTANCECFE_H_
#define _POINTVARIPLANEDISTANCECFE_H_

#include <Element.d/MpcElement.d/ConstraintFunction.d/PointVariPlaneDistanceConstraintFunction.h>
#include <Element.d/MpcElement.d/ConstraintFunctionElement.h>

class PointVariPlaneDistanceConstraintElement : public ConstraintFunctionElement<PointVariPlaneDistanceConstraintFunction>
{
  public:
    PointVariPlaneDistanceConstraintElement(int* _nn); 

  protected:
    void getConstants(CoordSet& cs, Eigen::Array<double,17,1>& sconst, Eigen::Array<int,1,1>& iconst, GeomState* = NULL);
};

#endif
