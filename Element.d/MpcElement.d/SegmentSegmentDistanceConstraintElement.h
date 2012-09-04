#ifndef _SEGMENTSEGMENTDISTANCECFE_H_
#define _SEGMENTSEGMENTDISTANCECFE_H_

#include <Element.d/MpcElement.d/ConstraintFunction.d/SegmentSegmentDistanceConstraintFunction.h>
#include <Element.d/MpcElement.d/ConstraintFunctionElement.h>

class SegmentSegmentDistanceConstraintElement : public ConstraintFunctionElement<SegmentSegmentDistanceConstraintFunction>
{
    double x1[3], x2[3]; // coordinates of the 2 points defining the line

  public:
    SegmentSegmentDistanceConstraintElement(int* _nn); 
    void setFrame(EFrame *);

  protected:
    void getConstants(CoordSet& cs, Eigen::Array<double,17,1>& sconst, Eigen::Array<int,1,1>& iconst, GeomState* = NULL);
};

#endif
