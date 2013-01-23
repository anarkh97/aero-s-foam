#ifndef _SPHERICALJOINT_H_
#define _SPHERICALJOINT_H_

#include <Element.d/Joint.d/BuildingBlocks.d/CommonPointConstraint.h>

class SphericalJoint : public CommonPointConstraint
{
  public:
    SphericalJoint(int*);
    int getTopNumber();
    PrioInfo examine(int sub, MultiFront*);
};

#endif
