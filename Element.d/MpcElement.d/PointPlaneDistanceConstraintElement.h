#ifndef _POINTPLANEDISTANCECFE_H_
#define _POINTPLANEDISTANCECFE_H_

#include <Element.d/MpcElement.d/ConstraintFunction.d/PointPlaneDistanceConstraintFunction.h>
#include <Element.d/MpcElement.d/ConstraintFunctionElement.h>

class PointPlaneDistanceConstraintElement : public ConstraintFunctionElement<PointPlaneDistanceConstraintFunction>
{
    double x1[3], x2[3], x3[3]; // coordinates of the 3 points defining the plane

  public:
    PointPlaneDistanceConstraintElement(int* _nn); 
    void setFrame(EFrame *);

  protected:
    void getConstants(CoordSet& cs, Eigen::Array<double,17,1>& sconst, Eigen::Array<int,1,1>& iconst, GeomState* = NULL);
};

#endif