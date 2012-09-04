#ifndef _POINTPOINTDISTANCECFE_H_
#define _POINTPOINTDISTANCECFE_H_

#include <Element.d/MpcElement.d/ConstraintFunction.d/PointPointDistanceConstraintFunction.h>
#include <Element.d/MpcElement.d/ConstraintFunctionElement.h>

class PointPointDistanceConstraintElement : public ConstraintFunctionElement<PointPointDistanceConstraintFunction>
{
    double x1[3]; // coordinates of the fixed point

  public:
    PointPointDistanceConstraintElement(int* _nn); 
    void setFrame(EFrame *);

  protected:
    void getConstants(CoordSet& cs, Eigen::Array<double,11,1>& sconst, Eigen::Array<int,1,1>& iconst, GeomState* = NULL);
};

#endif
