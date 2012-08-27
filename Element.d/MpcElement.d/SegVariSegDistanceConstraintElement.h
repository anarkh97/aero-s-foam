#ifndef _SEGVARISEGDISTANCECFE_H_
#define _SEGVARISEGDISTANCECFE_H_

#include <Element.d/MpcElement.d/ConstraintFunction.d/SegVariSegDistanceConstraintFunction.h>
#include <Element.d/MpcElement.d/ConstraintFunctionElement.h>

class SegVariSegDistanceConstraintElement : public ConstraintFunctionElement<SegVariSegDistanceConstraintFunction>
{

  public:
    SegVariSegDistanceConstraintElement(int* _nn); 

  protected:
    void getConstants(CoordSet& cs, Eigen::Array<double,17,1>& sconst, Eigen::Array<int,1,1>& iconst, GeomState* = NULL);
};

#endif
