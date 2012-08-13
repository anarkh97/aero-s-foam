#ifndef _PRISMATICJOINT_H_
#define _PRISMATICJOINT_H_

#include <Element.d/SuperElement.h>

// reference: Rigid Body Dynamics of Mechanisms: Theoretical basis, Volume 1
// Hubert Hahn, section 5.2.2.7
// constrains three rotational and two translational dofs

class PrismaticJoint : public SuperElement
{
  public:
    PrismaticJoint(int*);
    int getTopNumber();
    bool hasRot() { return true; }
    PrioInfo examine(int sub, MultiFront*);
};

#endif
