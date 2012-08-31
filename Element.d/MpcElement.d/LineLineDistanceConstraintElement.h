#ifndef _LINELINEDISTANCECFE_H_
#define _LINELINEDISTANCECFE_H_

#include <Element.d/MpcElement.d/ConstraintFunction.d/LineLineDistanceConstraintFunction.h>
#include <Element.d/MpcElement.d/ConstraintFunctionElement.h>

class LineLineDistanceConstraintElement : public ConstraintFunctionElement<LineLineDistanceConstraintFunction>
{
    double x1[3], x2[3]; // coordinates of the 2 points defining the line

  public:
    LineLineDistanceConstraintElement(int* _nn); 
    void setFrame(EFrame *);

  protected:
    void getConstants(CoordSet& cs, Eigen::Array<double,17,1>& sconst, Eigen::Array<int,1,1>& iconst, GeomState* = NULL);
};

#endif
