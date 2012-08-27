#ifndef _PLANARJOINT_H_
#define _PLANARJOINT_H_

#include <Element.d/SuperElement.h>

// constrains one translational and two rotational dofs

class PlanarJoint : public SuperElement
{
  public:
    PlanarJoint(int*);
    int getTopNumber();
    bool hasRot() { return true; }
    PrioInfo examine(int sub, MultiFront*);
};

#endif
