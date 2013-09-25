#ifndef _LINEVARILINEDISTANCECONSTRAINTELEMENT_H_
#define _LINEVARILINEDISTANCECONSTRAINTELEMENT_H_

#include <Element.d/Function.d/Constraint.d/LineVariLineDistanceConstraintFunction.h>
#include <Element.d/MpcElement.d/ConstraintFunctionElement.h>

class LineVariLineDistanceConstraintElement : public ConstraintFunctionElement<Simo::LineVariLineDistanceConstraintFunction>
{

  public:
    LineVariLineDistanceConstraintElement(int* _nn); 

  protected:
    void getConstants(CoordSet& cs, Eigen::Array<double,17,1>& sconst, Eigen::Array<int,1,1>& iconst, GeomState* = NULL);
};

#endif
