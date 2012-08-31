#ifndef _DISTANCECFE_H_
#define _DISTANCECFE_H_

#include <Element.d/MpcElement.d/ConstraintFunction.d/DistanceConstraintFunction.h>
#include <Element.d/MpcElement.d/ConstraintFunctionElement.h>

class DistanceConstraintElement : public ConstraintFunctionElement<DistanceConstraintFunction>
{
  public:
    DistanceConstraintElement(int* _nn, double f0, int type = 0); 

  protected:
    double f0;
    void getConstants(CoordSet& cs, Eigen::Array<double,4,1>& sconst, Eigen::Array<int,0,1>&, GeomState* = NULL);
};

#endif
