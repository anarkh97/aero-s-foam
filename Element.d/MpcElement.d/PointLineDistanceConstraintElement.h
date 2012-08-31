#ifndef _POINTLINEDISTANCECFE_H_
#define _POINTLINEDISTANCECFE_H_

#include <Element.d/MpcElement.d/ConstraintFunction.d/PointLineDistanceConstraintFunction.h>
#include <Element.d/MpcElement.d/ConstraintFunctionElement.h>

class PointLineDistanceConstraintElement : public ConstraintFunctionElement<PointLineDistanceConstraintFunction>
{
    double x1[3], x2[3]; // coordinates of the 2 points defining the line

  public:
    PointLineDistanceConstraintElement(int* _nn); 
    void setFrame(EFrame *);

  protected:
    void getConstants(CoordSet& cs, Eigen::Array<double,14,1>& sconst, Eigen::Array<int,1,1>& iconst, GeomState* = NULL);
};

#endif
