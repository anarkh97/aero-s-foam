#ifndef _DISTANCECONSTRAINTELEMENT_H_
#define _DISTANCECONSTRAINTELEMENT_H_

#include <Element.d/Function.d/Constraint.d/DistanceConstraintFunction.h>
#include <Element.d/MpcElement.d/ConstraintFunctionElement.h>

class DistanceConstraintElement : public ConstraintFunctionElement<Simo::DistanceConstraintFunction>
{
  public:
    DistanceConstraintElement(int* _nn, double f0, int type = 0); 
    double getVelocityConstraintRhs(GeomState*, GeomState&, CoordSet&, double);
    double getAccelerationConstraintRhs(GeomState*, GeomState&, CoordSet&, double);

  protected:
    double f0;
    void getConstants(CoordSet& cs, Eigen::Array<double,4,1>& sconst, Eigen::Array<int,0,1>&, GeomState* = NULL);
};

#endif
